#!/usr/bin/env ruby

require 'fileutils'
require 'tmpdir'

root = File.expand_path('../../../..', __dir__)
script = File.join(root, 'contrib/bin/medchem_transform_activity.rb')
fixture = __dir__

Dir.mktmpdir('medchem_transform_activity_test') do |tmpdir|
  output = File.join(tmpdir, 'output.txt')
  dbname = File.join(tmpdir, 'structures.bdb')

  cmd = [
    RbConfig.ruby,
    script,
    '-A', File.join(fixture, 'activity.txt'),
    '-d', dbname,
    '-o', output,
    '--buildsmidb', File.join(fixture, 'fake_buildsmidb.rb'),
    '--medchem-wizard', File.join(fixture, 'fake_medchem_wizard.rb'),
    '--in-database', File.join(fixture, 'fake_in_database.rb'),
    '--model-command', File.join(fixture, 'fake_model.rb'),
    File.join(fixture, 'input.smi')
  ]

  unless system(*cmd)
    warn "Command failed: #{cmd.join(' ')}"
    exit 1
  end

  expected = File.read(File.join(fixture, 'expected.txt'))
  got = File.read(output)
  if got != expected
    warn 'Output mismatch'
    warn 'Expected:'
    warn expected
    warn 'Got:'
    warn got
    exit 1
  end
end

puts 'medchem_transform_activity test passed'
