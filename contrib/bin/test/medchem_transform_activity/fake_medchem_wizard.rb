#!/usr/bin/env ruby

smiles = ARGV.fetch(-1)

File.foreach(smiles) do |line|
  tokens = line.split
  next if tokens.empty?

  _smiles, id = tokens
  case id
  when 'mol1'
    puts "CC mol1 add_carbon"
    puts "CO mol1 add_oxygen"
  when 'mol2'
    puts "CCC mol2 add_carbon"
    puts "CN mol2 add_nitrogen"
  end
end
