#!/usr/bin/env ruby
require 'open3'
exe, indir = ARGV
File.write('rectangular.txt', File.read(File.join(indir, 'rectangular.txt')))
File.write('ragged.txt', "a b\nc\n")
File.write('empty.txt', '')

def run(exe, args, success: true, input: '')
  out, err, status = Open3.capture3(exe, *args, stdin_data: input)
  abort "#{args.inspect}: status #{status.exitstatus}: #{err}" unless status.success? == success
  err
end

run(exe, ['rectangular.txt'])
run(exe, ['empty.txt'])
run(exe, ['-'], input: "a b\nc d\n")
run(exe, ['missing.txt'], success: false)
[['missing.txt', 'rectangular.txt'], ['rectangular.txt', 'missing.txt']].each do |args|
  err = run(exe, args, success: false)
  abort 'Valid file was not processed' unless err.include?("2 records read from 'rectangular.txt'")
end
# Keep shape diagnostics informational, including early termination with -b.
[[], ['-g'], ['-b'], ['-r']].each do |options|
  run(exe, [*options, 'ragged.txt'])
end
run(exe, ['ragged.txt', 'empty.txt', 'missing.txt'], success: false)
# The failure flag must not wrap back to zero with many missing inputs.
run(exe, Array.new(256, 'missing.txt'), success: false)
run(exe, ['-e', '9999', 'rectangular.txt'], success: false)
Dir.mkdir('unreadable_as_file')
run(exe, ['unreadable_as_file'], success: false)
run(exe, ['-C', 'nonexistent/counts.txt', 'rectangular.txt'], success: false)
puts 'tcount exit-status checks passed'
