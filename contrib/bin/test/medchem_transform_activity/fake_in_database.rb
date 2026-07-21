#!/usr/bin/env ruby

found = nil
notfound = nil
args = ARGV.dup

until args.empty?
  arg = args.shift
  found = args.shift if arg == '-F'
  notfound = args.shift if arg == '-U'
end

raise 'missing -F' unless found
raise 'missing -U' unless notfound

File.open(found, 'w') do |output|
  output.puts 'CC mol1 add_carbon mol2'
  output.puts 'CCC mol2 add_carbon mol3'
end

File.open(notfound, 'w') do |output|
  output.puts 'CO mol1 add_oxygen'
  output.puts 'CN mol2 add_nitrogen'
end

STDIN.read
