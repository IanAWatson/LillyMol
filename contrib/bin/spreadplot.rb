#!/usr/bin/env ruby

$stderr << "Once upon a time there was a tool to generate a plot of distances from spread.\n"
$stderr << "That tool does not need to exist, since making the plot is easy\n"
$stderr << "Run spread to get an output file\n"
$stderr << " gfp_spread file.gfp > file.spr\n"
$stderr << "Run nplotnn so you get the distance on each line and zero neighbours.\n";
$stderr << " nplotnn -c 4 -n 0 -s file.spr > file.spr.smi\n"
$stderr << "That file contains a column of distances - take that to your favourite\n"
$stderr << "plotting package\n"
