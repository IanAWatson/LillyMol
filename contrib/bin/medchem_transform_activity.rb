#!/usr/bin/env ruby

require 'csv'
require 'English'
require 'fileutils'
require 'open3'
require 'optparse'
require 'shellwords'
require 'tmpdir'

Stats = Struct.new(:n, :min, :max, :mean, :median, keyword_init: true)

class Values
  def initialize
    @values = []
  end

  def add(value)
    @values << value
  end

  def stats
    return Stats.new(n: 0, min: nil, max: nil, mean: nil, median: nil) if @values.empty?

    sorted = @values.sort
    nvalues = sorted.length
    median = if nvalues.odd?
               sorted[nvalues / 2]
             else
               0.5 * (sorted[nvalues / 2 - 1] + sorted[nvalues / 2])
             end

    Stats.new(
      n: nvalues,
      min: sorted.first,
      max: sorted.last,
      mean: @values.sum / nvalues,
      median:
    )
  end
end

def die(message)
  warn message
  exit 1
end

def quote_command(argv)
  argv.map(&:to_s).shelljoin
end

def run_command(argv, verbose:)
  warn quote_command(argv) if verbose
  system(*argv)
  die "Command failed: #{quote_command(argv)}" unless $CHILD_STATUS.success?
end

def run_pipeline(commands, verbose:)
  warn commands.map { |cmd| quote_command(cmd) }.join(' | ') if verbose

  wait_threads = Open3.pipeline_start(*commands)
  wait_threads.each_with_index do |thread, ndx|
    status = thread.value
    next if status.success?

    die "Pipeline command #{ndx + 1} failed: #{quote_command(commands[ndx])}"
  end
end

def each_non_blank_line(fname)
  return enum_for(__method__, fname) unless block_given?

  File.foreach(fname).with_index(1) do |line, line_number|
    line.chomp!
    next if line.empty?

    yield line, line_number
  end
end

def read_smiles_ids(fname)
  ids = {}
  each_non_blank_line(fname) do |line, line_number|
    tokens = line.split
    die "#{fname}:#{line_number}: expected 2 tokens, got #{tokens.size}" unless tokens.size == 2

    _smiles, id = tokens
    die "#{fname}:#{line_number}: duplicate id '#{id}'" if ids.key?(id)

    ids[id] = true
  end

  die "#{fname}: no molecules read" if ids.empty?
  ids
end

def read_csv_activity(fname)
  activity = {}

  CSV.foreach(fname, headers: true).with_index(2) do |row, line_number|
    die "#{fname}:#{line_number}: expected 2 columns, got #{row.fields.size}" unless row.fields.size == 2

    id = row[0].to_s
    die "#{fname}:#{line_number}: duplicate id '#{id}'" if activity.key?(id)

    activity[id] = Float(row[1])
  rescue ArgumentError
    die "#{fname}:#{line_number}: invalid activity '#{row[1]}'"
  end

  activity
end

def read_space_activity(fname)
  activity = {}
  saw_header = false

  each_non_blank_line(fname) do |line, line_number|
    unless saw_header
      saw_header = true
      next
    end

    tokens = line.split
    die "#{fname}:#{line_number}: expected 2 columns, got #{tokens.size}" unless tokens.size == 2

    id, value = tokens
    die "#{fname}:#{line_number}: duplicate id '#{id}'" if activity.key?(id)

    activity[id] = Float(value)
  rescue ArgumentError
    die "#{fname}:#{line_number}: invalid activity '#{value}'"
  end

  activity
end

def read_activity(fname, ids)
  activity = if File.extname(fname).casecmp?('.csv')
               read_csv_activity(fname)
             else
               read_space_activity(fname)
             end

  die "#{fname}: no activity values read" if activity.empty?

  ids.each_key do |id|
    die "#{fname}: missing activity for smiles id '#{id}'" unless activity.key?(id)
  end

  activity
end

def accumulate_found(fname, activity, observed)
  each_non_blank_line(fname) do |line, line_number|
    tokens = line.split
    die "#{fname}:#{line_number}: expected 4 tokens, got #{tokens.size}" unless tokens.size == 4

    _smiles, starting_id, transformation, found_id = tokens
    die "#{fname}:#{line_number}: unknown starting id '#{starting_id}'" unless activity.key?(starting_id)
    die "#{fname}:#{line_number}: unknown found id '#{found_id}'" unless activity.key?(found_id)

    observed[transformation].add(activity[found_id] - activity[starting_id])
  end
end

# Optional predicted-value input. This deliberately mirrors the not-found stream
# with one extra value appended by whatever model/prediction wrapper is used:
#   smiles starting_id transformation predicted_activity
def accumulate_predictions(fname, activity, predicted)
  each_non_blank_line(fname) do |line, line_number|
    tokens = line.split
    die "#{fname}:#{line_number}: expected 4 tokens, got #{tokens.size}" unless tokens.size == 4

    _smiles, starting_id, transformation, value = tokens
    die "#{fname}:#{line_number}: unknown starting id '#{starting_id}'" unless activity.key?(starting_id)

    predicted[transformation].add(Float(value) - activity[starting_id])
  rescue ArgumentError
    die "#{fname}:#{line_number}: invalid predicted activity '#{value}'"
  end
end


def command_from_template(command, input, output)
  argv = Shellwords.split(command)
  die 'Empty model command' if argv.empty?

  used_input = false
  used_output = false
  argv = argv.map do |token|
    updated = token.gsub('{input}') do
      used_input = true
      input
    end
    updated.gsub('{output}') do
      used_output = true
      output
    end
  end

  argv << input unless used_input
  argv << output unless used_output
  argv
end

def format_stat(value)
  return '.' if value.nil?

  format('%.6g', value)
end

def write_table(output, observed, predicted)
  output.puts %w[
    transformation
    observed_n observed_min observed_max observed_mean observed_median
    predicted_n predicted_min predicted_max predicted_mean predicted_median
  ].join(' ')

  transformations = (observed.keys + predicted.keys).uniq.sort
  transformations.each do |transformation|
    observed_stats = observed[transformation].stats
    predicted_stats = predicted[transformation].stats

    output.puts [
      transformation,
      observed_stats.n,
      format_stat(observed_stats.min),
      format_stat(observed_stats.max),
      format_stat(observed_stats.mean),
      format_stat(observed_stats.median),
      predicted_stats.n,
      format_stat(predicted_stats.min),
      format_stat(predicted_stats.max),
      format_stat(predicted_stats.mean),
      format_stat(predicted_stats.median)
    ].join(' ')
  end
end

def usage(parser)
  warn <<~USAGE
    Accumulate activity changes by medchem_wizard transformation.

    #{parser}

    The smiles input must contain exactly two tokens per non-blank line:
      smiles id

    The activity file must contain a header and two columns:
      id activity

    If the activity file name ends in .csv, it is parsed as CSV.
    Otherwise it is parsed as whitespace separated text.
  USAGE
  exit 1
end

options = {
  buildsmidb: 'buildsmidb_bdb',
  in_database: 'in_database_bdb',
  medchem_wizard: 'medchem_wizard.sh',
  output: '-',
  verbose: false,
  keep: false
}

parser = OptionParser.new do |opts|
  opts.banner = "Usage: #{File.basename($PROGRAM_NAME)} -A activity -d dbname [options] input.smi"

  opts.on('-A', '--activity FILE', 'Activity file') { |value| options[:activity] = value }
  opts.on('-d', '--database NAME', 'BerkeleyDB database name') { |value| options[:database] = value }
  opts.on('-o', '--output FILE', 'Output table, default stdout') { |value| options[:output] = value }
  opts.on('--found FILE', 'Found structures file') { |value| options[:found] = value }
  opts.on('--notfound FILE', 'Not-found structures file') { |value| options[:notfound] = value }
  opts.on('--predictions FILE', 'Existing predictions: smiles starting_id transformation predicted_activity') do |value|
    options[:predictions] = value
  end
  opts.on('--model-command COMMAND', 'Command that converts notfound.smi to predictions') do |value|
    options[:model_command] = value
  end
  opts.on('--buildsmidb EXE', 'buildsmidb executable, default buildsmidb_bdb') do |value|
    options[:buildsmidb] = value
  end
  opts.on('--in-database EXE', 'in_database executable, default in_database_bdb') do |value|
    options[:in_database] = value
  end
  opts.on('--medchem-wizard EXE', 'medchem_wizard executable, default medchem_wizard.sh') do |value|
    options[:medchem_wizard] = value
  end
  opts.on('--keep', 'Keep temporary found/notfound files') { options[:keep] = true }
  opts.on('-v', '--verbose', 'Verbose execution') { options[:verbose] = true }
end

begin
  parser.parse!
rescue OptionParser::ParseError => e
  die e.message
end

usage(parser) unless ARGV.size == 1
usage(parser) unless options[:activity] && options[:database]
die 'Specify only one of --predictions and --model-command' if options[:predictions] && options[:model_command]

smiles = ARGV.fetch(0)
ids = read_smiles_ids(smiles)
activity = read_activity(options.fetch(:activity), ids)

tmpdir = nil
unless options[:found] && options[:notfound] && (options[:predictions] || !options[:model_command])
  tmpdir = Dir.mktmpdir('medchem_transform_activity')
  options[:found] ||= File.join(tmpdir, 'found.smi')
  options[:notfound] ||= File.join(tmpdir, 'notfound.smi')
  options[:predictions] ||= File.join(tmpdir, 'predictions.txt') if options[:model_command]
end

begin
  run_command(
    [options.fetch(:buildsmidb), '-d', options.fetch(:database), '-g', 'all', '-l', '-c', smiles],
    verbose: options.fetch(:verbose)
  )

  medchem = [options.fetch(:medchem_wizard), '-W', 'space', smiles]
  lookup = [
    options.fetch(:in_database), '-d', options.fetch(:database), '-l', '-c', '-p',
    '-F', options.fetch(:found), '-U', options.fetch(:notfound), '-i', 'smi', '-'
  ]
  run_pipeline([medchem, lookup], verbose: options.fetch(:verbose))

  observed = Hash.new { |hash, key| hash[key] = Values.new }
  predicted = Hash.new { |hash, key| hash[key] = Values.new }

  accumulate_found(options.fetch(:found), activity, observed)
  if options[:model_command]
    model = command_from_template(options.fetch(:model_command), options.fetch(:notfound),
                                  options.fetch(:predictions))
    run_command(model, verbose: options.fetch(:verbose))
  end
  accumulate_predictions(options.fetch(:predictions), activity, predicted) if options[:predictions]

  if options.fetch(:output) == '-'
    write_table($stdout, observed, predicted)
  else
    File.open(options.fetch(:output), 'w') { |file| write_table(file, observed, predicted) }
  end
ensure
  FileUtils.remove_entry(tmpdir) if tmpdir && !options.fetch(:keep)
end
