#!/usr/bin/env ruby

input = ARGV.fetch(0)
output = ARGV.fetch(1)

predictions = {
  'add_oxygen' => 3.0,
  'add_nitrogen' => 1.5
}

File.open(output, 'w') do |out|
  File.foreach(input) do |line|
    tokens = line.split
    next if tokens.empty?

    smiles, starting_id, transformation = tokens
    prediction = predictions.fetch(transformation)
    out.puts [smiles, starting_id, transformation, prediction].join(' ')
  end
end
