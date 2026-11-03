#!/usr/bin/env ruby
require 'open3'
exe, indir = ARGV
input = File.join(indir, 'input.smi')
def run(exe, input, args, success: true)
  out, err, status = Open3.capture3(exe, *args, input)
  abort "#{args.inspect}: #{err}" unless status.success? == success
  abort "Sanitizer failure: #{err}" if err.match?(/AddressSanitizer|runtime error:/)
  out
end
base = ['-s', 'C-Cl', '-k', '2', '-m', '1', '-M', '20']
plain = run(exe, input, base)
abort 'Default behaviour added caps' if plain.include?('*')
capped = run(exe, input, [*base, '-B', 'chiraldummy'])
expected = {'enantiomer1' => 'Br[C@@H](F)*', 'enantiomer2' => 'Br[C@H](F)*'}
expected.each do |name, smiles|
  abort "Lost orientation for #{name}: #{capped}" unless capped.include?("#{smiles} #{name} FRAGID=")
end
# Multiple cuts at one centre each add a separate unlabelled cap.
multi_args = ['-B', 'chiraldummy', '-s', 'C-[Cl,F,Br]', '-k', '3', '-m', '1', '-M', '20']
[[], ['-B', 'recap'], ['-B', 'lostchiral']].each do |extra|
  multiple = run(exe, input, [*multi_args, *extra])
  fragments = multiple.lines.grep(/FRAGID=/).map { |line| line.split.first }
  abort "Missing three-cap fragment: #{multiple}" unless fragments.any? { |smi| smi.count('*') == 3 }
  unless extra.include?('recap')
    abort "Missing two-cap fragment: #{multiple}" unless fragments.any? { |smi| smi.count('*') == 2 }
  end
  abort 'Dummy caps unexpectedly isotopic' if fragments.any? { |smi| smi.match?(/\[\d+\*/) }
  if extra.include?('lostchiral')
    abort 'Invalid chirality retained' if fragments.any? { |smi| smi.count('*') >= 2 && smi.include?('@') }
  end
end
# Reconstructing complements must also add every cap at the parent centre.
['def', 'auto'].each do |mode|
  text = run(exe, input, [*multi_args, '-C', mode])
  abort 'Text complement lacks multiple caps' unless text.split.any? { |token| token.count('*') >= 2 }
end
proto = run(exe, input, [*multi_args, '-C', 'def', '-B', 'proto'])
abort 'Proto complement lacks multiple caps' unless proto.scan(/(?:smi|comp): "([^"]+)"/).flatten.any? { |smi| smi.count('*') >= 2 }
# A four-connected centre can become the only real atom in a capped fragment.
File.write('four.smi', "F[C@](Cl)(Br)I four\n")
four = run(exe, 'four.smi', ['-B', 'chiraldummy', '-s', 'C-[F,Cl,Br,I]', '-k', '4', '-m', '1', '-M', '20'])
abort 'Missing four-cap central fragment' unless four.lines.grep(/FRAGID=/).any? { |line| line.split.first.count('*') == 4 }
four_clean = run(exe, 'four.smi', ['-B', 'chiraldummy', '-B', 'lostchiral', '-s', 'C-[F,Cl,Br,I]', '-k', '4', '-m', '1', '-M', '20'])
abort 'Cleanup prevented four caps' unless four_clean.lines.grep(/FRAGID=/).any? { |line| line.split.first.count('*') == 4 && !line.split.first.include?('@') }
recap = run(exe, input, [*base, '-B', 'chiraldummy', '-B', 'recap'])
expected.each_value { |smiles| abort 'RECAP lost cap orientation' unless recap.include?(smiles) }
text = run(exe, input, [*base, '-B', 'chiraldummy', '-C', 'def'])
expected.each do |name, smiles|
  abort 'Text complement lacks cap' unless text.include?("Cl #{name} #{smiles} COMP #{name}")
end
proto = run(exe, input, [*base, '-B', 'chiraldummy', '-C', 'def', '-B', 'proto'])
expected.each_value { |smiles| abort 'Proto complement lacks cap' unless proto.include?(%Q(comp: "#{smiles}")) }
limited = run(exe, input, ['-B', 'chiraldummy', '-s', 'C-Cl', '-k', '2', '-m', '1', '-M', '3'])
abort 'Caps not counted in maximum fragment size' if limited.lines.grep(/FRAGID=/).any? { |line| line.include?('*') }
# A cut elsewhere in a capped fragment can still recurse, and null provenance
# must be handled safely during parent-index mapping.
File.write('chain.smi', "Cl[C@H](F)CCOC chain\n")
run(exe, 'chain.smi', ['-B', 'chiraldummy', '-s', 'C-Cl', '-s', 'C-O', '-k', '3', '-m', '1', '-M', '20'])
%w[BRB BRF BRS dearom].each do |option|
  run(exe, input, [*base, '-B', 'chiraldummy', '-B', option], success: false)
end
run(exe, input, [*base, '-B', 'chiraldummy', '-I', 'enva'], success: false)
run(exe, input, [*base, '-B', 'chiraldummy', '-C', 'fragatype', '-P', 'UST:Z'], success: false)
puts 'chiral dummy checks passed'
