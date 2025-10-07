#!/usr/bin/env ruby
# Detects numeric differences between two files.
# Each pair of lines is tokenised.
# Each token that can be parsed as a float is compared.

require "#{ENV['LILLYMOL_HOME']}/contrib/bin/lib/iwcmdline"

def usage
end

class Options 
  attr_reader :input_separator, :relative_error, :absolute_error, :skip_header
  def initialize(cl)
    @input_separator = if cl.option_present('i')
                         cl.value('i')
                       else
                         ' '
                       end
      
    @relative_error = if cl.option_present('r')
                        cl.value('r')
                      else
                        1.0e-05
                      end
    @absolute_error = if cl.option_present('a')
                        cl.value('a')
                      else
                        1.0e-05
                      end
    @skip_header = cl.option_present('h')

  end

  # Return true if `v1` and `v2` are numerically equivalent
  def same?(v1, v2)
    return true if (v1 - v2).abs <= @absolute_error

    v1 = v1.abs
    v2 = v2.abs
    if v1 < v2
      v1, v2 = v2, v1
    end
    # $stderr << "v1 #{v1} v2 #{v2}\n"
    return true if (v1 - v2) / v1 <= @relative_error

    return false
  end
end

# Tokenise the lines `line1` and `line2` and return true if all tokens
# are numerically equivalent.
def numerically_equivalent_lines(line1, line2, line_number, options)
  f1 = line1.chomp.split(options.input_separator)
  f2 = line2.chomp.split(options.input_separator)

  unless f1.size == f2.size
    $stderr << "Token count mismatch on line #{line_number}\n"
    return false
  end

  f1.each_with_index do |t1, ndx|
    t2 = f2[ndx]
    next if t1 == t2
    # $stderr << "Cmp #{t1} with #{t2}\n"
    next if /[A-D,F-Z,a-d,f-z]/.match(t1)
    next if /[A-D,F-Z,a-d,f-z]/.match(t2)
    next unless /[0-9]/.match(t1)
    next unless /[0-9]/.match(t2)
    begin
      v1 = Float(t1)
      v2 = Float(t2)
    rescue
      next
    end
    next if v1 == v2
    next if options.same?(v1, v2)
    $stderr << "Numeric mismatch #{t1} #{t2}\n"
    return false
  end

  return true
end

# Return true if all tokens in files `f1` and `f2` are the same.
def numerically_equivalent_files(f1, f2, options)
  f1.each_line.each_with_index do |line1, line_number|
    line2 = f2.readline
    next if line_number == 0 and options.skip_header
    return false unless numerically_equivalent_lines(line1, line2, line_number, options)
  end

  return true
end

# Return true if all tokens in files `fname1` and `fname2` are the same.
def numerically_equivalent_file_names(fname1, fname2, options)
  $stderr <<  fname1 << ' ' << fname2 << "\n"
  File.open(fname1, 'r') do |file1|
    File.open(fname2, 'r') do |file2|
      return false unless numerically_equivalent_files(file1, file2, options)
    end
  end

  return true
end

def main
  unless ARGV.size == 2
    $stderr << "Must specify exactly two files to be compared\n"
    usage
  end

  cl = IWCmdline.new("-v-i=s-r=float-a=float-h")

  options = Options.new(cl)

  rc = numerically_equivalent_file_names(ARGV[0], ARGV[1], options)
  if rc
    exit(0)
  else
    exit(1)
  end
end

main
