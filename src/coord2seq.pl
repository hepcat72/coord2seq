#!/usr/bin/perl -w
#PBS -l nodes=1:ppn=1

#getsubseq.pl
#Generated using perl_script_template.pl 1.33
#Robert W. Leach
#rwleach@ccr.buffalo.edu
#Created on 4/1/2008
#Center for Computational Research
#Copyright 2007

#To-do:
#1. Edit getsubseq.pl so that the coordinate file is what the output suffix is
#   appended to, unless I've got multiple coord files for one set of sequences,
#   but I think that's a more rare case, which can be circumvented by piping
#   them in.
#2. The code for getsubseq appears to enforce unique subseq IDs globally
#   across all the files.  I understand why it does that - so when aggregating
#   all sequences in an STDOUT file, they all have unique IDs, but I should add
#   an option to check that they are only unique per file.

#These variables (in main) are used by printVersion()
my $template_version_number = '1.34';
my $software_version_number = '1.22';

##
## Start Main
##

use strict;
use Getopt::Long;

#Declare & initialize variables.  Provide default values here.
my($outfile_suffix); #Not defined so a user can overwrite the input file
my @fasta_files                = ();
my @coord_files                = ();
my @input_coords               = ();
my $current_output_file        = '';
my $help                       = 0;
my $version                    = 0;
my $force                      = 0;
my $linear_flag                = 0;
my $start_col_num              = 1;
my $stop_col_num               = 2;
my @subseq_id_col_nums         = ();
my $dir_col_num                = 4;
my $parent_id_col_num          = 5;
my $comment_delimiter          = ';';
my @comment_col_nums           = ();
my $allow_relative_coords      = 0;
my $allow_zeroes               = 0;
my $no_wrap                    = 0;
my $pad_char                   = 'n';
my $skip_done                  = 0;
my $subseqid_delim             = ':';
my $mask_mode                  = 0;
my $mask                       = 'N';

#These variables (in main) are used by the following subroutines:
#verbose, error, warning, debug, printVersion, getCommand and usage
my $preserve_args = [@ARGV];  #Preserve the agruments for getCommand
my $verbose       = 0;
my $quiet         = 0;
my $DEBUG         = 0;

my $GetOptHash =
  {'i|fasta-file=s'        => sub {push(@fasta_files,        #REQUIRED unless
                                        sglob($_[1]))},      # <> is supplied
   '<>'                    => sub {if($_[0]=~/[\D ]/)        #REQUIRED unless
				     {push(@fasta_files,     # -i is supplied
					   sglob($_[0]))}
				   else{push(@input_coords,
					     sglob($_[0]))}},
   'f|coord-file=s'        => sub {push(@coord_files,        #REQUIRED unless
					sglob($_[1]))},      # -c is supplied
   'c|coords=s'            => sub {push(@input_coords,       #REQUIRED unless
					grep {/\d/}          # -f is supplied
					split(/\D+/,
					      $_[1]))},
   'b|start-col-num=s'     => \$start_col_num,               #REQUIRED [1] unls
                                                             # -c is supplied
   'e|stop-col-num=s'      => \$stop_col_num,                #REQUIRED [2] unls
                                                             # -c is supplied
   'a|allow-relative-coords!' => \$allow_relative_coords,    #OPTIONAL [Off]
   'use-zero-coords!'      => \$allow_zeroes,                #OPTIONAL [Off]
   'no-wrap!'              => \$no_wrap,                     #OPTIONAL [Off]
   'l|linear!'             => \$linear_flag,                 #OPTIONAL [Off]
   's|subseqid-col-num=s'  => sub {push(@subseq_id_col_nums, #OPTIONAL [3]
					grep {/\d/}
					split(/\D+/,
					      $_[1]))},
   'd|direction-col-num=s' => \$dir_col_num,            #OPTIONAL [4]
   'p|parentid-col-num=s'  => \$parent_id_col_num,      #OPTIONAL [5]
   'm|comment-col-nums=s'  => sub {push(@comment_col_nums,#REQUIRED unless -f
					grep {/\d/}       #  is supplied
					split(/\D+/,
					      $_[1]))},
   'subseqid-delimiter=s'  => \$subseqid_delim,         #OPTIONAL [:]
   'comment-delimiter=s'   => \$comment_delimiter,      #OPTIONAL [;]
   'mask-mode!'            => \$mask_mode,              #OPTIONAL [Off]
   'mask-character=s'      => \$mask,                   #OPTIONAL [N]
   'o|outfile-suffix=s'    => \$outfile_suffix,         #OPTIONAL [undef]
   'force!'                => \$force,                  #OPTIONAL [Off]
   'skip-done!'            => \$skip_done,              #OPTIONAL [Off]
   'v|verbose!'            => \$verbose,                #OPTIONAL [Off]
   'q|quiet!'              => \$quiet,                  #OPTIONAL [Off]
   'h|help!'               => \$help,                   #OPTIONAL [Off]
   'debug!'                => \$DEBUG,                  #OPTIONAL [Off]
   'version!'              => \$version,                #OPTIONAL [Off]
  };

#If there are no arguments and no files directed or piped in
if(scalar(@ARGV) == 0 && isStandardInputFromTerminal())
  {
    usage();
    exit(0);
  }

#Get the input options
GetOptions(%$GetOptHash);

#Print the debug mode (it checks the value of the DEBUG global variable)
debug("Debug mode on.");

#If the user has asked for help, call the help subroutine
if($help)
  {
    help();
    exit(0);
  }

#If the user has asked for the software version, print it
if($version)
  {
    printVersion();
    exit(0);
  }

#Check validity of verbosity options
if($verbose && $quiet)
  {
    $quiet = 0;
    error("You cannot supply verbose and quiet flags at the same time.");
    exit(1);
  }

#Put standard input into the input_files array if standard input has been redirected in
if(!isStandardInputFromTerminal())
  {
    push(@fasta_files,'-') unless(scalar(grep {$_ eq '-'} (@fasta_files,
							   @coord_files)));

    #Warn the user about the naming of the outfile when using STDIN
    if(defined($outfile_suffix))
      {warning("Input on STDIN detected along with an outfile suffix.  Your ",
	       "output file will be named STDIN$outfile_suffix")}
  }

#Make sure there is input
if(scalar(@fasta_files) == 0)
  {
    error("No fasta files detected.");
    usage(1);
    exit(2);
  }

#Make sure there is coordinate input
if(scalar(@input_coords) == 0 && scalar(@coord_files) == 0)
  {
    error("No coordinates supplied via either -f or -c.");
    usage(1);
    exit(3);
  }

#Make sure there is only one type of coordinate input
if(scalar(@input_coords) && scalar(@coord_files))
  {
    error("The -f and -c options cannot be used together.  Please supply one ",
	  "or the other, not both.");
    usage(1);
    exit(4);
  }

#Make sure that the number of files is correct
if(scalar(@coord_files) != 0                           && #Coord files supplied
   scalar(@fasta_files) != scalar(@coord_files)        && #not same # as fasta
   scalar(@fasta_files) != 1 && scalar(@coord_files) != 1)#neither is 1 file
  {
    error("The number of coordinate and fasta files do not match.  You may ",
	  "either have 1 fasta file and numerous coordinate files, 1 ",
	  "coordinate file and numerous fasta files, or the same number of ",
	  "coordinate and fasta files.");
    usage(1);
    exit(5);
  }

if($skip_done && $force)
  {
    error("--skip-done and --force are incompatible.");
    usage(1);
    quit(1);
  }

#Check to make sure previously generated output files won't be over-written
if(!$skip_done && !$force && defined($outfile_suffix))
  {
    my $existing_outfiles = [];
    foreach my $output_file (map {($_ eq '-' ? 'STDIN' : $_) . $outfile_suffix}
			     @fasta_files)
      {push(@$existing_outfiles,$output_file) if(-e $output_file)}

    if(scalar(@$existing_outfiles))
      {
	error("The output files: [@$existing_outfiles] already exist.  ",
	      "Use -f to force overwrite.  E.g.\n\t",
	      getCommand(1),' --force');
	exit(6);
      }
  }
elsif($skip_done)
  {
    my @newinfiles    = ();
    my @newcoordfiles = ();
    my $numinfiles    = scalar(@fasta_files);
    my $numcoordfiles = scalar(@coord_files);
    my $f_index       = 0;
    foreach my $infile (@fasta_files)
      {
	my $coordfile = ($numcoordfiles == $numinfiles ?
			 $coord_files[$f_index] :
			 ($numcoordfiles ? $coord_files[0] : ''));
	if(-e (($infile eq '-' ? 'STDIN' : $infile) . $outfile_suffix))
	  {verbose("[$infile] Already done.  Skipping.")}
	else
	  {
	    push(@newinfiles,$infile);
	    push(@newcoordfiles,$coordfile)
	      if($numcoordfiles == $numinfiles);
	  }
	$f_index++;
      }
    @fasta_files = @newinfiles;
    @coord_files = @newcoordfiles;
  }

#Just to make sure no undefs got in there due to user error, filter for valid
#numbers
@comment_col_nums = grep {/\d/} @comment_col_nums;
my $comment_check = {};
my @comment_col_num_errs =
  grep {my $result = ($_ < 1 && exists($comment_check->{$_}));
	$comment_check->{$_} = 1;$result} @comment_col_nums;
if(scalar(@comment_col_num_errs))
  {
    $comment_check = {};
    @comment_col_nums =
      grep {my $result = ($_ > 0 && !exists($comment_check->{$_}));
	    $comment_check->{$_} = 1;$result} @comment_col_nums;
    error("These invalid or duplicate comment column numbers (-m): ",
	  "[@comment_col_num_errs].  These numbers have been removed.");
  }

if(scalar(@subseq_id_col_nums) == 0)
  {
    @subseq_id_col_nums =
      grep {$_} ($parent_id_col_num,$start_col_num,$stop_col_num);
  }
@subseq_id_col_nums = grep {/\d/} @subseq_id_col_nums;
my $subseq_id_check = {};
my @subseq_id_col_num_errs =
  grep {my $result = ($_ < 1 && exists($subseq_id_check->{$_}));
	$subseq_id_check->{$_} = 1;$result} @subseq_id_col_nums;
if(scalar(@subseq_id_col_num_errs))
  {
    $subseq_id_check = {};
    @subseq_id_col_nums =
      grep {my $result = ($_ > 0 && !exists($subseq_id_check->{$_}));
	    $subseq_id_check->{$_} = 1;$result} @subseq_id_col_nums;
    error("These invalid or duplicate subseq ID column numbers (-s): ",
	  "[@subseq_id_col_num_errs].  These numbers have been removed.");
  }

#Check that users haven't specified conflicting column numbers and warn them if
#so.  (i.e. let them do it)
my $columns = {};
push(@{$columns->{$start_col_num}},    'start');
push(@{$columns->{$stop_col_num}},     'stop');
push(@{$columns->{$dir_col_num}},      'direction');
push(@{$columns->{$parent_id_col_num}},'parent sequence ID');
foreach my $colnumkey (keys(%$columns))
  {
    next if($colnumkey == 0);
    if(scalar(@{$columns->{$colnumkey}}) > 1)
      {warning("You have specified these values: [",
	       join(', ',@{$columns->{$colnumkey}}),
	       "] to all be in column: [$colnumkey] of your coordinate ",
	       "file.  This looks like a mistake, but I'm going to move ",
	       "forward as if that's what you intended.  Use the column ",
	       "number options (-b, -e, -s, -d, -p, & -m) to correct this.")}
  }

if(isStandardOutputToTerminal() && !defined($outfile_suffix))
  {verbose("NOTE: VerboseOverMe functionality has been altered to yield ",
	   "clean STDOUT output.")}

my $universal_coords = []; #[{START=>$start,STOP=>},...]
if(scalar(@input_coords))
  {
    #Make sure there are an even number of coords
    if(scalar(@input_coords) % 2)
      {
	warning("There should be an even number of alternating start and ",
		"stop coordinates supplied to -c.  If you don't have a stop ",
		"coordinate, supply a 0 to indicate the end of the sequence.");
	push(@input_coords,0);
      }

    #Break up the coords into starts and stops
    for(my $i = 0;$i < scalar(@input_coords);$i += 2)
      {push(@$universal_coords,
	    {START=>$input_coords[$i],STOP=>$input_coords[$i+1]})}

    #Sort the array from smallest to largest start coord
    @$universal_coords = sort {$a->{START} <=> $b->{START}} @$universal_coords;
  }

if($allow_zeroes && !$allow_relative_coords)
  {
    error("--use-zero-coords was supplied, but --allow-relative-coords was ",
	  "not supplied.  --use-zero-coords only works if ",
	  "--allow-relative-coords is turned on.");
    usage(1);
    exit(7);
  }

if($no_wrap && !$allow_relative_coords)
  {
    error("--no-wrap was supplied, but --allow-relative-coords was not ",
	  "supplied.  --no-wrap only works if --allow-relative-coords is ",
	  "turned on.");
    usage(1);
    exit(7);
  }

verbose("Run conditions: ",getCommand(1),"\n");

#If output is going to STDOUT instead of output files with different extensions
if(!defined($outfile_suffix))
  {verbose("[STDOUT] Opened for all output.")}

my $all_id = '__aLL__';
my($coords);  #{defline_or_defline_id => [[start,stop],[start,stop],...]}
              #define_or_defline_id = "__aLL__" for coords relating to all seqs
my $num_coord_files = scalar(@coord_files);
my $num_fasta_files = scalar(@fasta_files);
my $first_loop      = 1;
my $parent_id_check = {};

#Calculate minimum number of columns to expect (not including comments):
my $min_num_columns = (sort {$b <=> $a} ($start_col_num,
					 $stop_col_num,
					 @subseq_id_col_nums,
					 $dir_col_num,
					 $parent_id_col_num))[0];
my $max_num_columns = (sort {$b <=> $a} ($start_col_num,
					 $stop_col_num,
					 @subseq_id_col_nums,
					 $dir_col_num,
					 $parent_id_col_num,
					 @comment_col_nums))[0];

#For each input file
foreach my $fasta_file (@fasta_files)
  {
    #If an output file name suffix has been defined
    if(defined($outfile_suffix))
      {
	##
	## Open and select the next output file
	##

	#Set the current output file name
	$current_output_file = ($fasta_file eq '-' ? 'STDIN' : $fasta_file)
	  . $outfile_suffix;

	#Open the output file
	if(!open(OUTPUT,">$current_output_file"))
	  {
	    #Report an error and iterate if there was an error
	    error("Unable to open output file: [$current_output_file]\n$!");
	    next;
	  }
	else
	  {verboseOverMe("[$current_output_file] Opened output file.")}

	#Select the output file handle
	select(OUTPUT);
      }

    ##
    ## Retrieve the coordinate system for this fasta file
    ##

    #If coordinates were supplied on the command line
    if(scalar(@$universal_coords))
      {$coords = {$all_id=>$universal_coords}}
    else
      {
	#If there's only one fasta file, get all the coords from all files
	if($num_fasta_files == 1)
	  {$coords = getCoordHash($allow_zeroes,
				  $min_num_columns,
				  $max_num_columns,
				  @coord_files)
	     if($first_loop)}
	#Else if there's >1 coordinate file (implying equal to # fasta files)
	elsif($num_coord_files > 1)
	  {$coords = getCoordHash($allow_zeroes,
				  $min_num_columns,
				  $max_num_columns,
				  shift(@coord_files))}
	#Else if $coords hasn't been filled (implying only one coordinate file)
	elsif(!defined($coords))
	  {$coords = getCoordHash($allow_zeroes,
				  $min_num_columns,
				  $max_num_columns,
				  $coord_files[0])}

	#Record each of the parent keys from the coords hash
	foreach my $parent_key (keys(%$coords))
	  {$parent_id_check->{$parent_key} = 0
	     unless(exists($parent_id_check->{$parent_key}))}
      }
    $first_loop = 0;

    #Open the input file
    if(!open(INPUT,$fasta_file))
      {
	#Report an error and iterate if there was an error
	error("Unable to open input file: [$fasta_file]\n$!");
	next;
      }
    else
      {verboseOverMe("[",($fasta_file eq '-' ? 'STDIN' : $fasta_file),
		     "] Opened input file.")}

    my $line_num = 0;
    my $line     = '';
    my $defline  = '';
    my($seq);
    my $num_extracted = 0;
    my $redundancy_check = {};

    #For each line in the current input file
    while(getLine(*INPUT))
      {
	$line_num++;
	$line = $_;

	verboseOverMe("[",
		      ($fasta_file eq '-' ? 'STDIN' : $fasta_file),
		      '] Reading line: [',$line_num,'].');

	next if($line !~ /\S/ || $line =~ /^\s*#/);
	if($line =~ />/)
	  {
	    if($defline)
	      {
		verboseOverMe("DOING Sequence: $defline");

		my $solidseq = formatSequence($seq,0);
		chomp($solidseq);
		chomp($defline);

		$num_extracted += extractSeqs($defline,$solidseq,$coords,
					      $parent_id_check,$mask_mode,
					      $mask);
	      }
	    $defline = $line;

	    my $tmp_id = $defline;
	    $tmp_id =~ s/^\s*>\s*//;
	    $tmp_id =~ s/\s.*//;
	    if($tmp_id eq '')
	      {warning("No Defline ID on line: [$line_num] of file: ",
		       "[$fasta_file].  Universal coordinates will be used ",
		       "if some were supplied either via command line ",
		       "arguments of via coordinate file with no parent ",
		       "sequence ID.")}
	    elsif(exists($redundancy_check->{$tmp_id}))
	      {
		error("Two sequences found with the same ID: [$tmp_id] in ",
		      "fasta file: [$fasta_file].  The same pairs of ",
		      "coordinates will be used for each sequence.");
	      }

	    $redundancy_check->{$tmp_id} = 1;

	    undef($seq);
	  }
	elsif($line =~ /^([^\t]+?) *\t\s*(.*)/)
	  {
	    $defline = $1;
	    $seq     = $2;

	    print STDERR ("DOING Sequence: $defline\n") if($verbose);

	    my $solidseq = formatSequence($seq,0);
	    chomp($solidseq);
	    chomp($defline);

	    $num_extracted += extractSeqs($defline,$solidseq,$coords,
					  $parent_id_check,$mask_mode,$mask);

	    debug("Running tally of sequences extracted: [$num_extracted].");

	    undef($seq);
	  }
	else
	  {$seq .= $line}
      }

    #Handle the last sequence (if there were any sequences)
    if(defined($seq))
      {
	my $solidseq = formatSequence($seq,0);
	chomp($solidseq);
	chomp($defline);

	$num_extracted += extractSeqs($defline,$solidseq,$coords,
				      $parent_id_check,$mask_mode,$mask);
      }

    debug("Total sequences extracted: [$num_extracted].");

    close(INPUT);

    if($num_extracted > 0)
      {verbose("$num_extracted sequences extracted.") unless($mask_mode)}
    else
      {
	if(-z $fasta_file)
	  {warning("[",($fasta_file eq '-' ? 'STDIN' : $fasta_file),
		   "] is empty.")}
	else
	  {error("No sequences were extracted for file [",
		 ($fasta_file eq '-' ? 'STDIN' : $fasta_file),"].")
	     unless($mask_mode)}
      }

    verbose("[",
	    ($fasta_file eq '-' ? 'STDIN' : $fasta_file),
	    '] Input file done.  Time taken: [',
	    scalar(markTime()),
	    " Seconds].");

    #If an output file name suffix is set
    if(defined($outfile_suffix))
      {
	#Select standard out
	select(STDOUT);
	#Close the output file handle
	close(OUTPUT);

	verbose("[$current_output_file] Output file done.");
      }
  }

my @missing = grep {$parent_id_check->{$_} == 0} (keys(%$parent_id_check));
if(scalar(@missing))
  {warning("These source sequences from the fasta files: [",
           join(',',sort @missing),
           "] were not involved in any subsequence extractions.")}

#Report the number of errors, warnings, and debugs
verbose("Done.  EXIT STATUS: [",
	"ERRORS: ",
	($main::error_number ? $main::error_number : 0),
	" WARNINGS: ",
	($main::warning_number ? $main::warning_number : 0),
	($DEBUG ?
	 " DEBUGS: " . ($main::debug_number ? $main::debug_number : 0) : ''),
        " TIME: ",scalar(markTime(0)),"s]");
if($main::error_number || $main::warning_number)
  {verbose("Scroll up to inspect errors and warnings.")}

##
## End Main
##






























##
## Subroutines
##

sub extractSeqs
  {
    my $defline       = $_[0];        #Defline from fasta file
    my $seq           = $_[1];        #Sequence under defline from fasta file
    my $coords_hash   = $_[2];        #Coordinate file that accompanied the
                                      #fasta file or it could be coords from
                                      #the command line
    my $parent_check  = $_[3];
    my $mask_mode     = $_[4];
    my $mask          = $_[5];
    my $universal_id  = $all_id;      #global variable
    my $linearity     = $linear_flag; #global variable
    my $seq_length    = length($seq);
    my $num_extracted = 0;

    my $id = $defline;               #First non-white-space string
    $id =~ s/^\s*>\s*//;             #  (minus the '>')
    $id =~ s/\s.*$//;
    $id = $all_id if($id eq '');

    my $effective_id = '';

    #See if this sequence has specific coords to retrieve
    if(exists($coords_hash->{$id}))
      {$effective_id = $id}
    elsif(exists($coords_hash->{$defline}))
      {$effective_id = $defline}
    elsif(exists($coords_hash->{$universal_id}))
      {$effective_id = $universal_id}
    else
      {return($num_extracted)}

    foreach my $indiv_coord_hash (@{$coords_hash->{$effective_id}})
      {
	my $start_pad       = '';
	my $stop_pad        = '';
	my $effective_start = 1;

	#Catch the case where the coords are both off one end in the no-wrap
	#scenario
	if($allow_relative_coords && $no_wrap &&
	   exists($indiv_coord_hash->{START}) &&
	   exists($indiv_coord_hash->{STOP})  &&
	   (($indiv_coord_hash->{START} < 0 &&
	     $indiv_coord_hash->{STOP} < 0) ||
	    #Note, a 0 value means the end of the genome
	    (($indiv_coord_hash->{START} > $seq_length ||
	      $indiv_coord_hash->{START} == 0) &&
	     ($indiv_coord_hash->{STOP} > $seq_length ||
	      $indiv_coord_hash->{STOP} == 0))))
	  {
	    error("Both coordinates: [$indiv_coord_hash->{START},",
		  "$indiv_coord_hash->{STOP}] encountered for ",
		  "sequence: [$defline] are out of bounds on the same side ",
		  "and --no-wrap was supplied.  Skipping.");
	    next;
	  }

	if(exists($indiv_coord_hash->{START}))
	  {
	    if($indiv_coord_hash->{START} >= ($seq_length * 2))
	      {
		warning("Start coordinate [$indiv_coord_hash->{START}] is ",
			"greater than or equal to twice the size of the ",
			"sequence [$seq_length].  Setting to sequence ",
			"size.");
		$effective_start = $seq_length;
	      }
	    elsif(abs($indiv_coord_hash->{START}) >= ($seq_length * 2))
	      {
		warning("Start coordinate [$indiv_coord_hash->{START}] is ",
			"less than or equal to -1 * twice the size of the ",
			"sequence [$seq_length].  Setting to 1.");
		$effective_start = 1;
	      }
	    #Generate real existing coordinates assuming that negative numbers
	    #or numbers larger than the sequence assume a circular genome
	    #regardless of the linear flag
	    elsif($indiv_coord_hash->{START} < 0)
	      {
		if($allow_relative_coords)
		  {
		    if($no_wrap)
		      {
			$effective_start = 1;
			$start_pad =
			  $pad_char x abs($indiv_coord_hash->{START});
		      }
		    else
		      {
			$effective_start = $seq_length -
			  (abs($indiv_coord_hash->{START}) % $seq_length) + 1;

			warning("Negative start coordinate: [",
				$indiv_coord_hash->{START},"] encountered ",
				"for sequence: [$defline].  Assuming a ",
				"circular coordinate system and changing the ",
				"start to: [$effective_start].");
		      }
		  }
		else
		  {
		    if($linear_flag)
		      {
			warning("Negative start coordinate: [",
				$indiv_coord_hash->{START},"] encountered ",
				"for sequence: [$defline].  Setting to 1 ",
				"since the source sequence is linear.");
			$effective_start = 1;
		      }
		    else
		      {
			error("Negative start coordinate: [",
			      $indiv_coord_hash->{START},"] encountered for ",
			      "sequence: [$defline].  Skipping.  Use -a to ",
			      "allow negative coordinates.");
			next;
		      }
		  }
	      }
	    elsif($indiv_coord_hash->{START} > $seq_length)
	      {
		if($allow_relative_coords)
		  {
		    if($no_wrap)
		      {
			$effective_start = $seq_length;
			$start_pad =
			  $pad_char x ($indiv_coord_hash->{START} -
				       $seq_length);
		      }
		    else
		      {
			$effective_start = $indiv_coord_hash->{START} %
			  $seq_length;

			warning("Start coordinate: [",
				$indiv_coord_hash->{START},"] is larger than ",
				"the source sequence: [$seq_length] ",
				"encountered for sequence: [$defline].  ",
				"Assuming a circular coordinate system and ",
				"changing the start to: [$effective_start].");
		      }
		  }
		else
		  {
		    if($linear_flag)
		      {
			warning("Start coordinate: [",
				$indiv_coord_hash->{START},"] is larger than ",
				"the source sequence: [$seq_length bp] ",
				"encountered for sequence: [$defline].  ",
				"Setting to  [$seq_length bp] since the ",
				"source sequence is linear.");
			$effective_start = $seq_length;
		      }
		    else
		      {
			error("Start coordinate: [",
			      $indiv_coord_hash->{START},"] is larger than ",
			      "the source sequence: [$seq_length bp] ",
			      "encountered for sequence: [$defline].  ",
			      "Skipping.  Use -a to allow coordinates larger ",
			      "than the genome.");
			next;
		      }
		  }
	      }
	    elsif($indiv_coord_hash->{START} == 0)
	      {$effective_start = $seq_length}
	    else
	      {$effective_start = $indiv_coord_hash->{START}}
	  }

	my $effective_stop = $seq_length;
	if(exists($indiv_coord_hash->{STOP}))
	  {
	    if($indiv_coord_hash->{STOP} >= ($seq_length * 2))
	      {
		warning("Stop coordinate [$indiv_coord_hash->{STOP}] is ",
			"greater than or equal to twice the size of the ",
			"sequence [$seq_length].  Setting to sequence ",
			"size.");
		$effective_stop = $seq_length;
	      }
	    elsif(abs($indiv_coord_hash->{STOP}) >= ($seq_length * 2))
	      {
		warning("Stop coordinate [$indiv_coord_hash->{STOP}] is less ",
			"than or equal to -1 * twice the size of the ",
			"sequence [$seq_length].  Setting to 1.");
		$effective_stop = 1;
	      }
	    elsif($indiv_coord_hash->{STOP} < 0)
	      {
		if($allow_relative_coords)
		  {
		    if($no_wrap)
		      {
			$effective_stop = 1;
			$stop_pad =
			  $pad_char x abs($indiv_coord_hash->{STOP});
		      }
		    else
		      {
			$effective_stop = $seq_length -
			  (abs($indiv_coord_hash->{STOP}) % $seq_length) + 1;

			warning("Negative stop coordinate: [",
				$indiv_coord_hash->{STOP},"] encountered for ",
				"sequence: [$defline].  Assuming a circular ",
				"coordinate system and changing the stop to: ",
				"[$effective_stop].");
		      }
		  }
		else
		  {
		    if($linear_flag)
		      {
			warning("Negative stop coordinate: [",
				$indiv_coord_hash->{STOP},"] encountered for ",
				"sequence: [$defline].  Setting to 1 ",
				"since the source sequence is linear.");
			$effective_stop = 1;
		      }
		    else
		      {
			error("Negative stop coordinate: [",
			      $indiv_coord_hash->{STOP},"] encountered for ",
			      "sequence: [$defline].  Skipping.  Use -a to ",
			      "allow negative coordinates.");
			next;
		      }
		  }
	      }
	    elsif($indiv_coord_hash->{STOP} > $seq_length)
	      {
		if($allow_relative_coords)
		  {
		    if($no_wrap)
		      {
			$effective_stop = $seq_length;
			$stop_pad =
			  $pad_char x ($indiv_coord_hash->{STOP} -
				       $seq_length);
		      }
		    else
		      {
			$effective_stop = $indiv_coord_hash->{STOP} %
			  $seq_length;

			warning("Stop coordinate: ",
				"[$indiv_coord_hash->{STOP}] is larger than ",
				"the source sequence: [$seq_length bp] ",
				"encountered for sequence: [$defline].  ",
				"Assuming a circular coordinate system and ",
				"changing the stop to: [$effective_stop].");
		      }
		  }
		else
		  {
		    if($linear_flag)
		      {
			warning("Stop coordinate: [$indiv_coord_hash->{STOP}]",
				" is larger than the source sequence: ",
				"[$seq_length bp] encountered for sequence: ",
				"[$defline].  Setting to $seq_length since ",
				"the source sequence is linear.");
			$effective_stop = $seq_length;
		      }
		    else
		      {
			error("Stop coordinate: [$indiv_coord_hash->{STOP}] ",
			      "is larger than the source sequence: ",
			      "[$seq_length bp] encountered for sequence: ",
			      "[$defline].  Skipping.  Use -a to allow ",
			      "coordinates larger than the genome.");
			next;
		      }
		  }
	      }
	    elsif($indiv_coord_hash->{STOP} == 0)
	      {
		#The effective start will be the end of the sequence, set above
	      }
	    else
	      {$effective_stop = $indiv_coord_hash->{STOP}}
	  }

	my $effective_forward = 1;
	if(exists($indiv_coord_hash->{DIR}) &&
	   $indiv_coord_hash->{DIR} eq 'reverse')
	  {
	    $effective_forward = 0;
	    if($linear_flag && $effective_start < $effective_stop)
	      {
		warning("Start/Stop coordinates appear to be unordered.  ",
			"Since the chromosome is linear, the start/stop for ",
			"sequence [",
			($indiv_coord_hash->{ID} ne '' ?
			 $indiv_coord_hash->{ID} :
			 "$effective_id $indiv_coord_hash->{START}-" .
			 "$indiv_coord_hash->{STOP}"),
			"] will be switched.") unless($mask_mode);
		($effective_start,$effective_stop) = ($effective_stop,
						      $effective_start);
	      }
	  }
	elsif(exists($indiv_coord_hash->{DIR}) &&
	      $indiv_coord_hash->{DIR} eq 'forward')
	  {
	    $effective_forward = 1;
	    if($linear_flag && $effective_start > $effective_stop)
	      {
		warning("Start/Stop coordinates appear to be unordered.  ",
			"Since the chromosome is linear, the start/stop for ",
			"sequence [",
			($indiv_coord_hash->{ID} ne '' ?
			 $indiv_coord_hash->{ID} :
			 "$effective_id $indiv_coord_hash->{START}-" .
			 "$indiv_coord_hash->{STOP}"),
			"] will be ",
			"switched.") unless($mask_mode);
		($effective_start,$effective_stop) = ($effective_stop,
						      $effective_start);
	      }
	  }
	else
	  {$effective_forward = areCoordsForward($effective_start,
						 $effective_stop,
						 $seq_length,
						 $linearity)}

	debug("Subseqid: [",(defined($indiv_coord_hash->{ID}) &&
			     $indiv_coord_hash->{ID} ne '' ?
			     $indiv_coord_hash->{ID} :
			     "$effective_id $indiv_coord_hash->{START}-" .
			     "$indiv_coord_hash->{STOP}"),"] is [",
	      ($effective_forward ? 'forward' : 'reverse'),
	      "] from [$effective_start] to [$effective_stop].");

	my $effective_id = $defline;
	if(exists($indiv_coord_hash->{ID}) && $indiv_coord_hash->{ID} ne '')
	  {$effective_id = ">$indiv_coord_hash->{ID}"}

	my $subseq = getCoordSeq($seq,
				 $effective_start,
				 $effective_stop,
				 $effective_forward,
				 1,
				 $linear_flag,
				 $mask_mode,
				 $mask);

	#If we're in mask mode, the subseq returned was the whole sequence with
	#the sequence in the coordinates masked, so save it in $seq so it can
	#be printed later
	if($mask_mode)
	  {$seq = $subseq}
	else
	  {
	    #If relative coordinates were used and no wrapping is allowed, pad
	    #the sequence with the pad charcter - as many as the coordinates go
	    #off the ends
	    if($allow_relative_coords && $no_wrap)
	      {$subseq = $start_pad . $subseq . $stop_pad}

	    print((defined($effective_id) && $effective_id ne '' ?
		   $effective_id : '>Sequence'),
		  " $effective_start..$effective_stop",
		  ($effective_forward ? '' : 'c'),
		  (!defined($indiv_coord_hash->{COMMENT}) ||
		   $indiv_coord_hash->{COMMENT} eq '' ? '' :
		   " $indiv_coord_hash->{COMMENT}"),
		  "\n",
		  formatSequence($subseq,60));
	  }

	my $id = $defline;
	$id =~ s/^\s*>\s*(\S+).*/$1/;

	if(exists($parent_check->{$id}))
	  {$parent_check->{$id}++}

	$num_extracted++;
      }

    if($mask_mode)
      {print("$defline\n",formatSequence($seq,60),"\n")}

    return($num_extracted);
  }


#This sub reads a coordinate file (as described by --help) and returns this
#hash structure:
#{$sourceseqid=>[{START=>$start,STOP=>$stop,DIR=>$dir,ID=>$subseqid},...],...}
sub getCoordHash
  {
    my $coord_hash    = {};                  #This is what's returned
    my $universal_id  = $all_id;             #global variable
    my $start_col     = $start_col_num;      #global variable
    my $stop_col      = $stop_col_num;       #global variable
    my @subseqid_cols = @subseq_id_col_nums; #global variable
    my $dir_col       = $dir_col_num;        #global variable
    my $parentid_col  = $parent_id_col_num;  #global variable
    my @comment_cols  = @comment_col_nums;   #global variable
    my $delimiter     = $comment_delimiter;  #global variable

    if(scalar(@_) < 1)
      {
	error("Not enough parameters: [",scalar(@_),"].  Must be at least 1.");
	return($coord_hash);
      }

    my $allow_zeroes = shift(@_);

    if(-e $allow_zeroes)
      {
	unshift(@_,$min_num_columns);
	warning("Invalid allow_zeroes flag sent in: [$allow_zeroes].  If a ",
		"file by that name exists, it will be assumed to be a ",
		"coordinate file.  Setting allow_zeroes to 0.");
	$allow_zeroes = 0;
      }

    my $min_num_columns   = shift(@_);

    if($min_num_columns !~ /^\d+$/)
      {
	if(-e $min_num_columns)
	  {unshift(@_,$min_num_columns)}
	warning("Invalid minimum number of columns sent in: ",
		"[$min_num_columns].  Setting to 1.");
	$min_num_columns = 1;
      }

    my $max_num_columns   = shift(@_);

    if($max_num_columns !~ /^\d+$/)
      {
	if(-e $max_num_columns)
	  {unshift(@_,$max_num_columns)}
	warning("Invalid maximum number of columns sent in: ",
		"[$max_num_columns].  Disabling maximum.");
	$max_num_columns = 0;
      }

    my $coord_files       = [@_];
    my $ids_found         = 0;         #If both of these turn out to be 1,
    my $empty_ids_found   = 0;         #there's a problem
    my($start,$stop,$subseqid,$dir,$strand,$sourceseqid,$comment);

    my $subseq_id_check = {};

    #For each input file
    foreach my $coord_file (@$coord_files)
      {
	#Open the coord file
	if(!open(COORD,$coord_file))
	  {
	    #Report an error and iterate if there was an error
	    error("Unable to open coord file: [$coord_file]\n$!");
	    next;
	  }
	else
	  {verboseOverMe("[",
			 ($coord_file eq '-' ? 'STDIN' : $coord_file),
			 "] Opened coord file.")}

	#Keep track of coord file's line number
	my $line_num = 0;

	#For each line in the current coord file
	while(getLine(*COORD))
	  {
	    $line_num++;
	    verboseOverMe("[",
			  ($coord_file eq '-' ? 'STDIN' : $coord_file),
			  "] Reading line $line_num.");

	    next if(/^\s*$/ || /^\s*#/);

	    $subseqid    = '';
	    $sourceseqid = '';
	    $start       = '';
	    $stop        = '';
	    $strand      = '';
	    $dir         = '';
	    $comment     = '';

	    chomp;
	    s/ +$//;
	    s/^ +//;

	    #Split on tabs (with possible errant white space surrounding them
	    #and on ".." in case the start and stop coordinates are written
	    #like this: 1..1000
	    my @cols = split(/(?: *\t *|(?<=\d)\.\.(?=[\d\-\+]))/,$_);

	    #If the number of columns is too few$seq_length
	    if(scalar(@cols) < $min_num_columns)
	      {
		#Add a split on 2 or more spaces
		@cols = split(/(?: *\t *|  +|(?<=\d)\.\.(?=[\d\-\+]))/,$_);

		if(scalar(@cols) > $max_num_columns)
		  {error("Line [$line_num]: [$_] is not parsing to the ",
			 "correct number of columns.  Please check it to ",
			 "make sure tab-characters are being used between ",
			 "all columns.")}
	      }

	    @cols = split(/\s+/,$_) unless(scalar(@cols) > 1);

	    $start       = (scalar(@cols) >= $start_col &&
			    defined($cols[$start_col - 1]) ?
			    $cols[$start_col - 1] : '');
	    $stop        = (scalar(@cols) >= $stop_col &&
			    defined($cols[$stop_col - 1]) ?
			    $cols[$stop_col - 1] : '');
	    $subseqid    = join($subseqid_delim,
				map {my $subseqid_col = $_;
				     ($subseqid_col > 0 &&
				      scalar(@cols) >= $subseqid_col &&
				      defined($cols[$subseqid_col - 1]) ?
				      $cols[$subseqid_col - 1] : '')}
				@subseqid_cols);
	    $strand      = ($dir_col > 0 &&
			    scalar(@cols) >= $dir_col &&
			    defined($cols[$dir_col - 1]) ?
			    $cols[$dir_col - 1] : '');
	    $sourceseqid = ($parentid_col > 0 &&
			    scalar(@cols) >= $parentid_col &&
			    defined($cols[$parentid_col - 1]) ?
			    $cols[$parentid_col - 1] : $universal_id);
	    if(scalar(@comment_cols))
	       {$comment = join($delimiter,
				map {($comment_cols[$_] > 0 &&
				      scalar(@cols) >= $comment_cols[$_] &&
				      defined($cols[$comment_cols[$_] - 1]) ?
				      $cols[$comment_cols[$_] - 1] : '')}
				(0..$#comment_cols))}

	    if($start =~ /[^\d\-\+\s]/ || $stop =~ /[^\-\+\d\s]/)
	      {
		error("Non-numeric characters were found in the ",
		      "start[$start] or stop[$stop] column(s) on line: ",
		      "[$line_num] of coordinate file: [$coord_file]: [$_].  ",
		      "Only positive integers are allowed.  Assuming this is ",
		      "a comment line & skipping.");
		next;
	      }

	    #Make sure there's no bogus characters
	    ($start,$stop) = grep {/\d/} split(/[^\d\-]+/,"$start $stop");

	    if($allow_zeroes)
	      {
		if($start < 1)
		  {$start--}
		if($stop < 1)
		  {$stop--}
	      }

	    #If the negative sign isn't first or the coord is 0
	    if($start =~ /\d\-/ || #$start == 0 ||
	       $stop  =~ /\d\-/) #|| $stop  == 0)
	      {
		error("An invalid start [$start] or stop [$stop] was found ",
		      "on line: [$line_num] of file: [$coord_file]: [$_].  ",
		      "The coordinate system is assumed to start at 1.  I'm ",
		      "going to assume this is a comment line & skip it.  ",
		      "See --use-zero-coords in the usage output if your ",
		      "coordinate file uses 0 to represent the last base in ",
		      "a circular genome or the position before the first ",
		      "base (when --no-wrap is used).");
		next;
	      }

	    #If there was no start, supply a default of 1
	    if($start eq '')
	      {$start = 1}

	    if($stop eq '')
	      {$stop = 0}

	    if($strand =~ /c|-|rev/i && $strand =~ /\+|for/i)
	      {
		error("Ambiguous strand: [$strand] on line: [$line_num].  It ",
		      "matches both the forward pattern: [+|for] and the ",
		      "reverse pattern: [c|-|rev].  The direction will be ",
		      "determined dynamically.");
		$dir = ($start < $stop || $stop == 0 ? 'forward' : 'reverse');
	      }
	    elsif($strand =~ /c|-|rev/i)
	      {$dir = 'reverse'}
	    elsif($strand =~ /\+|for|\A[123]\Z/i)
	      {$dir = 'forward'}
	    elsif($strand eq '')
	      {$dir = 'dynamic'}
	    else
	      {
		error("Unrecognized strand string: [$strand] on line: ",
		      "[$line_num].  Please use: ['forward', 'reverse', 'c', ",
		      "'+', or '-'].  See --help for more info.  The ",
		      "direction will be determined dynamically.");
		$dir = 'dynamic';
	      }

	    debug("Subseqid: [$subseqid] is [$dir].");

	    #Sort the start and stop
#	    ($start,$stop) = sort {$a <=> $b} ($start,$stop)
#	      unless($stop == 0);

	    #Trim the '>' if they put it there
	    $subseqid =~ s/^\s*>\s*//;
	    #Remove everything including and beyond the first space
	    $subseqid =~ s/ +.*$//;

	    #If no ID was parsed
	    if($subseqid !~ /\S/)
	      {
		$subseqid = '';
		$empty_ids_found = 1;
	      }
	    else
	      {
		$ids_found = 1;

		$subseq_id_check->{$subseqid}++;
	      }

	    #Populate the coord hash
	    push(@{$coord_hash->{$sourceseqid}},
		 {START=>$start,STOP=>$stop,DIR=>$dir,ID=>$subseqid,
		  COMMENT=>$comment});
	  }

	close(COORD);

	if(scalar(grep {$subseq_id_check->{$_} > 1} keys(%$subseq_id_check)))
	  {
	    warning("Found multiple instances of subsequence IDs: [",
		    join(',',grep {$subseq_id_check->{$_} > 1}
			 keys(%$subseq_id_check)),
		    "] in coordinate file: [$coord_file].  ",
		    "Output deflines will not have unique IDs.  ",
		    "(Note, this can be a side-effect if you do not ",
		    "submit the same number of Fasta and coordinate ",
		    "files.)");
	  }

	verbose("[",
		($coord_file eq '-' ? 'STDIN' : $coord_file),
		'] Coord file done.  Time taken: [',
		scalar(markTime()),
		" Seconds].");
      }

    return($coord_hash);
  }

#Copied from seq-lib.pl-new on 4/1/2008
sub getCoordSeq
  {
    #1. Read in the parameters.
    my $sequence          = $_[0];
    my $start             = $_[1];
    my $stop              = $_[2];
    my $forward_flag      = $_[3];
    my $force_smaller_seq = $_[4];  #Use when $start and $stop are not ordered
    my $linear_flag       = $_[5];
    my $mask_mode         = $_[6];
    my $mask              = ($_[7] ? $_[7] : 'N');
    my $sub_sequence;

    #2. Error check the parameters and fill in default values if optional
    #   parameters are not supplied and get rid of white space characters in
    #   the sequence.
    $sequence =~ s/\s+//g;        #White spaces are allowed
    $sequence =~ s/<[^>]*>//g;    #HTML tags are allowed
    if(length($sequence) < $start || length($sequence) < $stop)
      {
        error("The sequence length: [",
	      length($sequence),
	      "] is smaller than your coordinate(s): [$start,$stop]!");
        return '';
      }
    if(!defined($forward_flag))
      {$forward_flag = 1}
    #end if(!defined($forward_flag))

    #Swap start and stop, based on directionality to generate the smaller of
    #two possible sequences in a circular genome (either spanning the origin of
    #the parent sequence or not)
    if(!$linear_flag && defined($force_smaller_seq) && $force_smaller_seq)
      {
	my $slen = length($sequence);
	if($forward_flag)
	  {
	    if(($start > $stop &&
		abs($stop - $start) < (($slen - $start + 1) + $stop)) ||
	       ($stop > $start &&
		abs($stop - $start) > (($slen - $stop + 1) + $start)))
	      {($start,$stop) = ($stop,$start)}
	  }
	else
	  {
	    #If the stop is greater than the start or if the stop is less than
	    #the start and the length would be ridiculously big
	    if(($stop > $start &&
		abs($stop - $start) < (($slen - $stop + 1) + $start)) ||
	       ($stop < $start &&
		abs($stop - $start) > (($slen - $start + 1) + $stop)))
	      {($start,$stop) = ($stop,$start)}
	  }
      }

    debug("Getting [",($forward_flag ? 'forward' : 'reverse'),
	  "] sequence from start: [$start] to stop: [$stop].");

    if(defined($start) && ($start !~ /^\d+$/ || $start < 1))
      {
        error("The start: [$start] was either invalid or not supplied!");
        return '';
      }
    #end if(defined($start) && ($start !~ /^\d+$/ || $start < 1))
    elsif(!defined($start) && $forward_flag)
      {$start = 1}
    #end elsif($forward_flag)
    elsif(!defined($start))
      {$start = length($sequence)}
    #end else

    if(defined($stop) &&
       (($stop !~ /^\d+$/) || ($stop < 1) || ($stop > length($sequence))))
      {
        error("The stop: [$stop] is invalid!");
        return '';
      }
    #end if(defined($stop) &&
    #       (($stop !~ /^\d+$/) || ($stop < 1) || ($stop > length($sequence))))
    elsif(!defined($stop) && $forward_flag)
      {$stop = length($sequence)}
    #end elsif($forward_flag)
    elsif(!defined($stop))
      {$stop = 1}
    #end else

    #4. If forward_flag is true and start is greater than stop
    if($forward_flag && $start > $stop)
      {
	if($mask_mode)
	  {
	    substr($sequence,
		   $start-1,
		   length($sequence)-$start-1,
		   $mask x (length($sequence)-$start-1));
	    substr($sequence,0,$stop,$mask x $stop);
	  }
	else
	  {
	    #4.1 sub_sequence is set from the start to the end of the sequence
	    $sub_sequence  = substr($sequence,$start-1);
	    #4.2 The start of sequence to the stop is appended to sub_sequence
	    $sub_sequence .= substr($sequence,0,$stop);
	  }
      }
    #end if($forward_flag && $start > $stop)
    #5. Else if forward flag is false and start is less than stop
    elsif(!$forward_flag && $start < $stop)
      {
	if($mask_mode)
	  {
	    substr($sequence,
		   $stop-1,
		   length($sequence)-$stop-1,
		   $mask x (length($sequence)-$stop-1));
	    substr($sequence,0,$start,$mask x $start);
	  }
	else
	  {
	    #5.1 sub_sequence is set from the stop to the end of the sequence
	    $sub_sequence  = substr($sequence,$stop-1);
	    #5.2 The start of sequence to the start is appended to sub_sequence
	    $sub_sequence .= substr($sequence,0,$start);
	  }
      }
    #end elsif(!$forward_flag && $start < $stop)
    #6. Else if start is less than stop
    elsif($start <= $stop)
      {
	if($mask_mode)
	  {
	    substr($sequence,$start-1,($stop - $start + 1),
		   $mask x ($stop - $start + 1));
	  }
	else
	  {
	    #6.1 sub_sequence is set from the start to the stop
	    $sub_sequence = substr($sequence,$start-1,($stop - $start + 1));
	  }
      }
    #end elsif($start <= $stop)
    #7. Else
    else
      {
	if($mask_mode)
	  {
	    substr($sequence,$stop-1,($start - $stop + 1),
		   $mask x ($start - $stop + 1));
	  }
	else
	  {
	    #7.1 sub_sequence is set from the stop to the start
	    $sub_sequence = substr($sequence,$stop-1,($start - $stop + 1));
	  }
      }
    #end else

    #8. If forward_flag is false
    if(!$forward_flag && !$mask_mode)
      {
        #8.1 Return the reverse compliment of sub_sequence
        return reverseComplement($sub_sequence);
      }
    #end if(!$forward_flag)
    #9. Return the sub_sequence
    return($mask_mode ? $sequence : $sub_sequence);
  }

#Copied from seq-lib.pl-new on 4/1/2008
sub areCoordsForward
  {
    #1. Read in the parameters.
    my $start     = $_[0];
    my $stop      = $_[1];
    my $molsize   = $_[2];  #Optional: defaults to infinity
                            #          (i.e. treats as linear)
    my $linearity = $_[3];  #Optional: default: circular if molsize is supplied
                            #                   otherwise linear
                            #If the string is not 'linear' or 'circular', it's
                            #treated as a flag where non-zero is linear and
                            #otherwise circular

    my $linear_flag = 0;
    if($linearity =~ /linear/i)
      {$linear_flag = 1}
    elsif($linearity !~ /circ/i && $linearity)
      {$linear_flag = 1}
    if($start !~ /^\d+$/ || $stop !~ /^\d+$/ ||
       (defined($molsize) && $molsize !~ /^\d+$/))
      {
        error("Invalid or unsupplied start: [$start] or stop: [$stop] or ",
              "invalid molecule size: [$molsize]");
      }
    #end if($start !~ /^\d+$/ || stop !~ /^\d+$/ ||
    #       (defined(molsize) && molsize !~ /^\d+$/))
    #2. If molsize is not defined OR the absolute value of start minus stop is
    #   less than molsize/2
    if(!defined($molsize)                      ||
       (defined($linear_flag) && $linear_flag) ||
       abs($start - $stop) < $molsize/2)
      {
	debug("Start less than stop?");
        #2.1 Return the result of this expression: start <= stop
        return $start <= $stop;
      }
    debug("Start greater than stop?");
    #end if(!defined(molsize) || abs($start - stop) < molsize/2)
    #3. Return the result of this expression: start >= stop

    #Kludge - need to change this logic all-together! - 2014/04/15
    if($start > $stop && $molsize > 1000000)
      {
	my $strand = $start >= $stop ? 'forward' : 'reverse';
	warning("Since we are in circular genome mode, no strand was ",
		"provided, the distance between the start and stop is more ",
		"than half the size of the parent sequence, and the parent ",
		"sequence is big enough to be a chromosome (bacterial), ",
		"we're guessing that start [$start] and stop [$stop] refer ",
		"to a sequence on the [$strand] strand.  Provide -l or ",
		"strand information to avoid guessing.");
      }

    return($molsize > 1000000 ? ($start >= $stop) : ($start <= $stop));
  }


#Copied from seq-lib.pl on 9/9/04 so as to be independent -Rob
sub formatSequence
  {
    #1. Read in the parameters.
    my $sequence          = $_[0];
    my $chars_per_line    = $_[1];
    my $coords_left_flag  = $_[2];
    my $coords_right_flag = $_[3];
    my $start_coord       = $_[4];
    my $coords_asc_flag   = $_[5];
    my $coord_upr_bound   = $_[6];
    my $uppercase_flag    = $_[7];
    my $print_flag        = $_[8];
    my $nucleotide_flag   = $_[9];

    my($formatted_sequence,
       $sub_string,
       $sub_sequence,
       $coord,
       $max_num_coord_digits,
       $line_size_left,
       $lead_spaces,
       $line);
    my $coord_separator = '  ';
    my $tmp_sequence = $sequence;
    $tmp_sequence =~ s/\s+//g;
    $tmp_sequence =~ s/<[^>]*>//g;
    my $seq_len = length($tmp_sequence);

    #2. Error check the parameters and set default values if unsupplied.
    my $default_chars_per_line    = ''; #Infinity
    my $default_coords_left_flag  = 0;
    my $default_coords_right_flag = 0;
    my $default_start_coord       = (!defined($coords_asc_flag) ||
				     $coords_asc_flag ? 1 : $seq_len);
    my $default_coords_asc_flag   = 1;
    my $default_coord_upr_bound   = undef();  #infinity (going past 1 produces
    my $default_uppercase_flag    = undef();  #          negative numbers)
    my $default_print_flag        = 0;

    if(!defined($chars_per_line) || $chars_per_line !~ /^\d+$/)
      {
        if($chars_per_line !~ /^\d+$/ && $chars_per_line =~ /./)
	  {warning("Invalid chars_per_line: [$chars_per_line] - using ",
		   "default: [$default_chars_per_line].")}
        #end if(chars_per_line !~ /^\d+$/)
	$chars_per_line = $default_chars_per_line;
      }
    elsif(!$chars_per_line)
      {$chars_per_line = ''}
    #end if(!defined($chars_per_line) || $chars_per_line !~ /^\d+$/)
    if(!defined($coords_left_flag))
      {$coords_left_flag = $default_coords_left_flag}
    #end if(!defined(coords_left_flag))
    if(!defined($coords_right_flag))
      {$coords_right_flag = $default_coords_right_flag}
    #end if(!defined(coords_right_flag))
    if(!defined($start_coord) || $start_coord !~ /^\-?\d+$/)
      {
        if(defined($start_coord) &&
	   $start_coord !~ /^\d+$/ && $start_coord =~ /./ &&
           ($coords_left_flag || $coords_right_flag))
          {warning("Invalid start_coord: [$start_coord] - using default: ",
		   "[$default_start_coord].")}
        #end if($start_coord !~ /^\d+$/)
        $start_coord = $default_start_coord;
      }
    #end if(!defined($start_coord) || $start_coord !~ /^\d+$/)
    if(!defined($coords_asc_flag))
      {$coords_asc_flag = $default_coords_asc_flag}
    #end if(!defined(coords_right_flag))
    if(defined($coord_upr_bound) && $coord_upr_bound !~ /^\d+$/)
      {undef($coord_upr_bound)}
    if(!defined($print_flag))
      {$print_flag = $default_print_flag}
    #end if(!defined($print_flag))

    if(defined($coord_upr_bound) && $start_coord < 1)
      {$start_coord = $coord_upr_bound + $start_coord}
    elsif($start_coord < 1)
      {$start_coord--}
    elsif(defined($coord_upr_bound) && $start_coord > $coord_upr_bound)
      {$start_coord -= $coord_upr_bound}

    #3. Initialize the variables used for formatting.  (See the DATASTRUCTURES
    #   section.)
    if($coords_asc_flag)
      {
        if(defined($coord_upr_bound) &&
           ($seq_len + $start_coord) > $coord_upr_bound)
          {$max_num_coord_digits = length($coord_upr_bound)}
        else
          {$max_num_coord_digits = length($seq_len + $start_coord - 1)}

        $coord = $start_coord - 1;
      }
    else
      {
        if(defined($coord_upr_bound) && ($start_coord - $seq_len + 1) < 1)
          {$max_num_coord_digits = length($coord_upr_bound)}
        elsif(!defined($coord_upr_bound) &&
              length($start_coord - $seq_len - 1) > length($start_coord))
          {$max_num_coord_digits = length($start_coord - $seq_len - 1)}
        else
          {$max_num_coord_digits = length($start_coord)}

        $coord = $start_coord + 1;
      }
    $line_size_left = $chars_per_line;
    $lead_spaces    = $max_num_coord_digits - length($start_coord);

    #5. Add the first coordinate with spacing if coords_left_flag is true.
    $line = ' ' x $lead_spaces . $start_coord . $coord_separator
      if($coords_left_flag);

    #6. Foreach sub_string in the sequence where sub_string is either a
    #   sub_sequence or an HTML tag.
    foreach $sub_string (split(/(?=<)|(?<=>)/,$sequence))
      {
        #6.1 If the substring is an HTML tag
        if($sub_string =~ /^</)
          #6.1.1 Add it to the current line of the formatted_sequence
          {$line .= $sub_string}
        #end if(sub_string =~ /^</)
        #6.2 Else
        else
          {
            $sub_string =~ s/\s+//g;

	    if($nucleotide_flag)
	      {
		my(@errors);
		(@errors) = ($sub_string =~ /([^ATGCBDHVRYKMSWNX])/ig);
		$sub_string =~ s/([^ATGCBDHVRYKMSWNX])//ig;
		if(scalar(@errors))
		  {warning('[',
			   scalar(@errors),
			   "] bad nucleotide characters were filtered out of ",
			   "your sequence: [",
			   join('',@errors),
			   "].")}
	      }

            #6.2.1 If the sequence is to be uppercased
            if(defined($uppercase_flag) && $uppercase_flag)
              #6.2.1.1 Uppercase the sub-string
              {$sub_string = uc($sub_string)}
            #end if(defined($uppercase_flag) && $uppercase_flag)
            #6.2.2 Else if the sequence is to be lowercased
            elsif(defined($uppercase_flag) && !$uppercase_flag)
              #6.2.2.1 Lowercase the sub-string
              {$sub_string = lc($sub_string)}
            #end elsif(defined($uppercase_flag) && !$uppercase_flag)

            #6.2.3 While we can grab enough sequence to fill the rest of a line
            while($sub_string =~ /(.{1,$line_size_left})/g)
              {
                $sub_sequence = $1;
                #6.2.3.1 Add the grabbed sequence to the current line of the
                #        formatted sequence
                $line .= $sub_sequence;
                #6.2.3.2 Increment the current coord by the amount of sequence
                #        grabbed
                my $prev_coord = $coord;
                if($coords_asc_flag)
                  {
                    $coord += length($sub_sequence);
                    if(defined($coord_upr_bound)      &&
                       $prev_coord <= $coord_upr_bound &&
                       $coord > $coord_upr_bound)
                      {$coord -= $coord_upr_bound}
                  }
                else
                  {
                    $coord -= length($sub_sequence);
                    if(defined($coord_upr_bound) &&
                       $prev_coord >= 1 && $coord < 1)
                      {$coord = $coord_upr_bound + $coord - 1}
                    elsif($prev_coord >= 1 && $coord < 1)
                      {$coord--}
                  }
#print STDERR "max_num_coord_diogits: [$max_num_coord_digits] coord: [$coord] line_size_left: [$line_size_left] sub_sequence: [$sub_sequence]\n";
                #6.2.3.3 If the length of the current sequence grabbed
                #        completes a line
                if(($line_size_left eq '' && length($sub_sequence) == 0) ||
		   ($line_size_left ne '' &&
		    length($sub_sequence) == $line_size_left))
                  {
                    $lead_spaces = $max_num_coord_digits - length($coord);
                    #6.2.3.3.1 Conditionally add coordinates based on the
                    #          coords flags
                    $line .= $coord_separator . ' ' x $lead_spaces . $coord
                      if($coords_right_flag);

                    #6.2.3.3.2 Add a hard return to the current line of the
                    #          formatted sequence
                    $line .= "\n";

                    #6.2.3.3.3 Add the current line to the formatted_sequence
                    $formatted_sequence .= $line;
                    #6.2.3.3.4 Print the current line if the print_flag is true
                    print $line if($print_flag);

                    #6.2.3.3.5 Start the next line
                    $lead_spaces = $max_num_coord_digits - length($coord+1);
                    $line = '';
                    $line = ' ' x $lead_spaces
                          . ($coords_asc_flag ? ($coord+1) : ($coord-1))
                          . $coord_separator
                      if($coords_left_flag);

                    #6.2.3.3.6 Reset the line_size_left (length of remaining
                    #          sequence per line) to chars_per_line
                    $line_size_left = $chars_per_line;
                  }
                #end if(length($sub_sequence) == $line_size_left)
                #6.2.3.4 Else
                else
                  #6.2.3.4.1 Decrement line_size_left (length of remaining
                  #          sequence per line) by the amount of sequence
                  #          grabbed
                  {
		    $line_size_left = 0 if($line_size_left eq '');
		    $line_size_left -= length($sub_sequence);
		  }
                #end 6.2.3.4 Else
              }
            #end while($sub_string =~ /(.{1,$line_size_left})/g)
          }
        #end 6.2 Else
      }
    #end foreach $sub_string (split(/(?=<)|(?<=>)/,$sequence))
    #7. Add the last coodinate with enough leadin white-space to be lined up
    #   with the rest coordinates if the coords_right_flag is true
    $lead_spaces = $max_num_coord_digits - length($coord);
    $line .= ' ' x $line_size_left . $coord_separator . ' ' x $lead_spaces
          . $coord
      if($coords_right_flag && $line_size_left != $chars_per_line);
    $line =~ s/^\s*\d+$coord_separator\s*$// if($coords_left_flag);

    #8. Add the ending PRE tag to the last line of the formatted sequence
    $line =~ s/\n*$/\n/s;

    #9. Add the last line to the formatted_sequence
    $formatted_sequence .= $line;
    #10. Print the last line if the print_flag is true
    print $line if($print_flag);

    if($coord < 1 && ($coords_left_flag || $coords_right_flag))
      {warning("The sequence straddles the origin.  Coordinates are ",
	       "inaccurate.")}

    #11. Return the formatted_sequence
    return $formatted_sequence;
  }

#copied from seq-lib.pl on 11/8/2005
sub reverseComplement
  {
    #1. Read in the sequence parameter.
    my $sequence = $_[0];
    my @errors;
    if(@errors =
       ($sequence =~ /([^ATGCBVDHRYKMSWNatgcbvdhrykmswn\.\-\s\r])/isg))
      {warning("Bad character(s) found: ['",join("','",@errors),"'].")}
    #end if(@errors = ($sequence =~ /([^ATGCBVDHRYKM\s\r])/isg))
    #2. Translate the new_sequence.
    $sequence =~ tr/ATGCBVDHRYKMatgcbvdhrykm/TACGVBHDYRMKtacgvbhdyrmk/;
    #3. Reverse the new_sequence.
    $sequence = reverse($sequence);
    return $sequence;
  }



##
## This subroutine prints a description of the script and it's input and output
## files.
##
sub help
  {
    my $script = $0;
    my $lmd = localtime((stat($script))[9]);
    $script =~ s/^.*\/([^\/]+)$/$1/;

    #Print a description of this program
    print << "end_print";

$script
Copyright 2007
Robert W. Leach
Created on DATE HERE
Last Modified on $lmd
Center for Computational Research
701 Ellicott Street
Buffalo, NY 14203
rwleach\@ccr.buffalo.edu

* WHAT IS THIS: This script takes fasta files and coordinate files (or
                coordinates supplied on the command line) and extracts the
                sequences from the fasta file corresponding to the supplied
                coordinates (inclusive and indexed from a start coordinate of
                1).

* FASTA FILE FORMAT: Standard fasta format: sequence lines preceded by a
                     defline which starts with a '>'.

* COORDINATE FILE FORMAT: White-space-delimited (e.g. tabs) file with the
                          columns start, stop, optional subsequence ID,
                          optional strand, and optional source sequence ID.
                          Spaces are not allowed other than to separate
                          columns.  The source sequence ID refers to the first
                          non-whitespace characters on the fasta defline of the
                          input sequence file.  If no sequence ID is provided,
                          the sequence is extracted from every supplied FASTA
                          entry.  Acceptable strand column values are: (for
                          forward strand:) '+' or 'forward' (for reverse
                          complement strand:) '-', 'reverse', or 'c'.  If no
                          subsequence ID is supplied, the defline of the source
                          sequence is used followed by the subsequence
                          coordinates (delimited by a space).  If no strand is
                          supplied the shorter sequence will be returned.  This
                          assumes a circular genome.  To avoid sequences
                          spanning the origin, supply the --linear flag.  If
                          the start or stop is not supplied, 1 and the end of
                          the sequence are used respectively.  Example:

                          1	300	Orf1	+	chrom1
                          301	499	Orf2	-	chrom1
                          907	634	Orf3		chrom1

                          In the above example, the strand for Orf3 is
                          determined dynamically as the reverse complement
                          strand because the start is larger than the stop.

* OUTPUT FORMAT: Fasta files with deflines formatted like this:

                   >Subsequence_ID <start>..<stop><c> <comment1;comment2>

                 Or like this (if you do not provide subsequence IDs):

                   >original_defline <start>..<stop><c> <comment1;comment2>

                 Note that the original defline of the source sequence is used
                 when no sebsequence ID is provided.

                 Comments are optional (-m) and the semicolon delimiter can be
                 changed using the --comment-delimiter option.

                 An input fasta file that has a defline like this:

                   >contig00227 length=213451 numreads=53105

                 will look like this when no comment columns are provided and a
                 coordinate entry asks for the reverse complement sequence from
                 coordinate 213361 to 212762:

                   >contig00227 length=213451 numreads=53105 213361..212762c


end_print

    return(0);
  }

##
## This subroutine prints a usage statement in long or short form depending on
## whether "no descriptions" is true.
##
sub usage
  {
    my $no_descriptions = $_[0];

    my $script = $0;
    $script =~ s/^.*\/([^\/]+)$/$1/;

    #Grab the first version of each option from the global GetOptHash
    my $options = '[' .
      join('] [',
	   grep {$_ ne '-i'}        #Remove REQUIRED params
	   map {my $key=$_;         #Save the key
		$key=~s/\|.*//;     #Remove other versions
		$key=~s/(\!|=.)$//; #Remove trailing getopt stuff
		$key = (length($key) > 1 ? '--' : '-') . $key;} #Add dashes
	   grep {$_ ne '<>'}        #Remove the no-flag parameters
	   keys(%$GetOptHash)) .
	     ']';

    print << "end_print";
USAGE: $script -i "input file(s)" $options
       $script $options < input_file
end_print

    if($no_descriptions)
      {print("Execute $script with no options to see a description of the ",
             "available parameters.\n")}
    else
      {
        print << 'end_print';

     -i|--fasta-file*     REQUIRED Space-separated fasta file(s inside quotes).
                                   *No flag required.  Standard input via
                                   redirection is acceptable.  Perl glob
                                   characters (e.g. '*') are acceptable inside
                                   quotes.  See --help for a description of
                                   fasta format.  Note, file names are not
                                   allowed to consist of only numbers because
                                   it will be interpretted as a command-line
                                   supplied coordinate.
     -f|--coord-file      REQUIRED Do not supply if -c has been supplied.
                                   Space-separated coordinate file(s inside
                                   quotes).  Perl glob characters (e.g. '*')
                                   are acceptable inside quotes.  See --help
                                   for a description of the coordinate file
                                   format.
     -c|--coords*         REQUIRED Do not supply if -f has been supplied.
                                   *No flag required.  Even number of (space-
                                   separated if no flag is used) coordinates.
                                   If the flag is used, the coordinates may be
                                   separated by non-numeric characters.  Note,
                                   negative coordinates are not allowed.
     -l|--linear          OPTIONAL [circular] Use this flag to treat submitted
                                   sequences as linear sequences instead of
                                   circular sequences (meaning no sequences
                                   will "span the origin").  Default behavior
                                   is for bacterial chromosomes and plasmids.
     -a|--allow-relative- OPTIONAL [Off] Allow coordinates that are larger than
        coords                     the source sequence (i.e. chromosome) and
                                   coordinates that are less than 1.  These
                                   coordinates will be recalculated assuming a
                                   circular genome.
     --use-zero-coords    OPTIONAL [Off] Makes the first coordinate of a
                                   sequence start at zero.  Must be supplied
                                   with -a.
     --no-wrap            OPTIONAL [Off] Prevents both start and stop
                                   coordinates from being recalculated when -a
                                   is supplied.  This will cause an error if,
                                   for example, both the start and stop are
                                   negative or if they are both larger than the
                                   genome.
     -b|--start-col-num   OPTIONAL [1] The column number in which to find the
                                   start coordinate from the input coordinate
                                   file(s).  Note that a start and stop
                                   separated by ".." are considered separate
                                   columns.  Columns are assumed to be
                                   delimited by tabs and instances of "..".
     -e|--stop-col-num    OPTIONAL [2] The column number in which to find the
                                   stop coordinate from the input coordinate
                                   file(s).  Note that a start and stop
                                   separated by ".." are considered separate
                                   columns.  Columns are assumed to be
                                   delimited by tabs and instances of "..".
     -s|--subseqid-col-   OPTIONAL [5,1,2*] The column number(s) in which to
        num                        find the sub-sequence ID from the input
                                   coordinate file(s).  If multiple subseqid
                                   columns are specified, the
                                   --subseqid-delimiter will delimit the
                                   portions of the ID, concatenated in the
                                   column order provided.  Note, this column is
                                   optional.  Note that columns are assumed to
                                   be delimited by tabs and instances of "..",
                                   so count start..stop as 2 columns when
                                   supplying this column number.
                                   * The default is the values of the -p, -b, &
                                   -e options supplied, in that order, which
                                   each have a default of 5, 1, & 2
                                   respectively.
     -d|--direction-col-  OPTIONAL [4] The column number in which to find the
        num                        direction (+/-,forward/reverse,c) from the
                                   input coordinatefile(s).  Note, this column
                                   is optional.  Default is the shorter
                                   sequence.  Also note that if you use 'c' to
                                   indicate reverse, you must explicitly
                                   specify 'forward' for all the rest.  Any
                                   missing direction will revert to default
                                   shorter sequence behavior.  If your
                                   coordinate file does not have a direction
                                   column, set this to 0 (or if there is no
                                   column 4, you can leave it as the default).
                                   Note that columns are assumed to be
                                   delimited by tabs and instances of "..", so
                                   count start..stop as 2 columns when
                                   supplying this column number.
     -p|--parentid-col-   OPTIONAL [5] The column number in which to find the
        num                        sub-sequence ID from the input coordinate
                                   file(s).  Note, this column is optional.  If
                                   your coordinate file does not have parent
                                   sequence IDs, set this to 0 (or if there is
                                   no column 5, you can leave it as the
                                   default).  Note that columns are assumed to
                                   be delimited by tabs and instances of "..",
                                   so count start..stop as 2 columns when
                                   supplying this column number.
     -m|--comment-col-    OPTIONAL [none] A non-number (e.g. space) delimited
        nums                       list of column numbers (inside quotes).
                                   Columns specified here will be appended to
                                   the defline of ">ID start..stop" with a
                                   delimiting space.  If multiple comment
                                   columns are specified, the
                                   --comment-delimiter will separate comment
                                   columns only.  Note that columns are assumed
                                   to be delimited by tabs and instances of
                                   "..", so count start..stop as 2 columns when
                                   supplying this column number.
     --subseqid-delimiter OPTIONAL [:] The delimiter to use between values
                                   found in the subseqid columns (-s) when
                                   appending them to the defline.  Note that if
                                   the column is empty, the delimiters will
                                   still be printed.
     --comment-delimiter  OPTIONAL [;] The delimiter to use between values
                                   found in the comment columns (-m) when
                                   appending them to the defline.  Note that if
                                   the column is empty, the delimiters will
                                   still be printed.
     --mask-mode          OPTIONAL [Off] Causes the sequences represented in
                                   the coordinates to be masked with the
                                   --mask-character instead of extracting the
                                   sub-sequences.
     --mask-character     OPTIONAL [N] If --mask-mode is supplied, this is the
                                   masking character that subsequences will be
                                   replaced with.
     -o|--outfile-suffix  OPTIONAL [nothing] This suffix is added to the input
                                   file names to use as output files.
                                   Redirecting a file into this script will
                                   result in the output file name to be "STDIN"
                                   with your suffix appended.
     --force              OPTIONAL [Off] Force overwrite of existing output
                                   files.  Only used when the -o option is
                                   supplied.  Cannot be used with skip-done.
     --skip-done          OPTIONAL Skip any input files that already have an
                                   output file.  Cannot be used with --force.
     -v|--verbose         OPTIONAL [Off] Verbose mode.  Cannot be used with the
                                   quiet flag.
     -q|--quiet           OPTIONAL [Off] Quiet mode.  Turns off warnings and
                                   errors.  Cannot be used with the verbose
                                   flag.
     -h|--help            OPTIONAL [Off] Help.  Use this option to see an
                                   explanation of the script and its input and
                                   output files.
     --version            OPTIONAL [Off] Print software version number.  If
                                   verbose mode is on, it also prints the
                                   template version used to standard error.
     --debug              OPTIONAL [Off] Debug mode.

end_print
      }

    return(0);
  }


##
## Subroutine that prints formatted verbose messages.  Specifying a 1 as the
## first argument prints the message in overwrite mode (meaning subsequence
## verbose, error, warning, or debug messages will overwrite the message
## printed here.  However, specifying a hard return as the first character will
## override the status of the last line printed and keep it.  Global variables
## keep track of print length so that previous lines can be cleanly
## overwritten.
##
sub verbose
  {
    return(0) unless($verbose);

    #Read in the first argument and determine whether it's part of the message
    #or a value for the overwrite flag
    my $overwrite_flag = $_[0];

    #If a flag was supplied as the first parameter (indicated by a 0 or 1 and
    #more than 1 parameter sent in)
    if(scalar(@_) > 1 && ($overwrite_flag eq '0' || $overwrite_flag eq '1'))
      {shift(@_)}
    else
      {$overwrite_flag = 0}

    #Ignore the overwrite flag if STDOUT will be mixed in
    $overwrite_flag = 0 if(isStandardOutputToTerminal());

    #Read in the message
    my $verbose_message = join('',@_);

    $overwrite_flag = 1 if(!$overwrite_flag && $verbose_message =~ /\r/);

    #Initialize globals if not done already
    $main::last_verbose_size  = 0 if(!defined($main::last_verbose_size));
    $main::last_verbose_state = 0 if(!defined($main::last_verbose_state));
    $main::verbose_warning    = 0 if(!defined($main::verbose_warning));

    #Determine the message length
    my($verbose_length);
    if($overwrite_flag)
      {
	$verbose_message =~ s/\r$//;
	if(!$main::verbose_warning && $verbose_message =~ /\n|\t/)
	  {
	    warning("Hard returns and tabs cause overwrite mode to not work ",
		    "properly.");
	    $main::verbose_warning = 1;
	  }
      }
    else
      {chomp($verbose_message)}

    if(!$overwrite_flag)
      {$verbose_length = 0}
    elsif($verbose_message =~ /\n([^\n]*)$/)
      {$verbose_length = length($1)}
    else
      {$verbose_length = length($verbose_message)}

    #Overwrite the previous verbose message by appending spaces just before the
    #first hard return in the verbose message IF THE VERBOSE MESSAGE DOESN'T
    #BEGIN WITH A HARD RETURN.  However note that the length stored as the
    #last_verbose_size is the length of the last line printed in this message.
    if($verbose_message =~ /^([^\n]*)/ && $main::last_verbose_state &&
       $verbose_message !~ /^\n/)
      {
	my $append = ' ' x ($main::last_verbose_size - length($1));
	unless($verbose_message =~ s/\n/$append\n/)
	  {$verbose_message .= $append}
      }

    #If you don't want to overwrite the last verbose message in a series of
    #overwritten verbose messages, you can begin your verbose message with a
    #hard return.  This tells verbose() to not overwrite the last line that was
    #printed in overwrite mode.

    #Print the message to standard error
    print STDERR ($verbose_message,
		  ($overwrite_flag ? "\r" : "\n"));

    #Record the state
    $main::last_verbose_size  = $verbose_length;
    $main::last_verbose_state = $overwrite_flag;

    #Return success
    return(0);
  }

sub verboseOverMe
  {verbose(1,@_)}

##
## Subroutine that prints errors with a leading program identifier containing a
## trace route back to main to see where all the subroutine calls were from,
## the line number of each call, an error number, and the name of the script
## which generated the error (in case scripts are called via a system call).
##
sub error
  {
    return(0) if($quiet);

    #Gather and concatenate the error message and split on hard returns
    my @error_message = split("\n",join('',@_));
    pop(@error_message) if($error_message[-1] !~ /\S/);

    $main::error_number++;

    my $script = $0;
    $script =~ s/^.*\/([^\/]+)$/$1/;

    #Assign the values from the calling subroutines/main
    my @caller_info = caller(0);
    my $line_num = $caller_info[2];
    my $caller_string = '';
    my $stack_level = 1;
    while(@caller_info = caller($stack_level))
      {
	my $calling_sub = $caller_info[3];
	$calling_sub =~ s/^.*?::(.+)$/$1/ if(defined($calling_sub));
	$calling_sub = (defined($calling_sub) ? $calling_sub : 'MAIN');
	$caller_string .= "$calling_sub(LINE$line_num):"
	  if(defined($line_num));
	$line_num = $caller_info[2];
	$stack_level++;
      }
    $caller_string .= "MAIN(LINE$line_num):";

    my $leader_string = "ERROR$main::error_number:$script:$caller_string ";

    #Figure out the length of the first line of the error
    my $error_length = length(($error_message[0] =~ /\S/ ?
			       $leader_string : '') .
			      $error_message[0]);

    #Put location information at the beginning of each line of the message
    foreach my $line (@error_message)
      {print STDERR (($line =~ /\S/ ? $leader_string : ''),
		     $line,
		     ($verbose &&
		      defined($main::last_verbose_state) &&
		      $main::last_verbose_state ?
		      ' ' x ($main::last_verbose_size - $error_length) : ''),
		     "\n")}

    #Reset the verbose states if verbose is true
    if($verbose)
      {
	$main::last_verbose_size = 0;
	$main::last_verbose_state = 0;
      }

    #Return success
    return(0);
  }


##
## Subroutine that prints warnings with a leader string containing a warning
## number
##
sub warning
  {
    return(0) if($quiet);

    $main::warning_number++;

    #Gather and concatenate the warning message and split on hard returns
    my @warning_message = split("\n",join('',@_));
    pop(@warning_message) if($warning_message[-1] !~ /\S/);

    my $leader_string = "WARNING$main::warning_number: ";

    #Figure out the length of the first line of the error
    my $warning_length = length(($warning_message[0] =~ /\S/ ?
				 $leader_string : '') .
				$warning_message[0]);

    #Put leader string at the beginning of each line of the message
    foreach my $line (@warning_message)
      {print STDERR (($line =~ /\S/ ? $leader_string : ''),
		     $line,
		     ($verbose &&
		      defined($main::last_verbose_state) &&
		      $main::last_verbose_state ?
		      ' ' x ($main::last_verbose_size - $warning_length) : ''),
		     "\n")}

    #Reset the verbose states if verbose is true
    if($verbose)
      {
	$main::last_verbose_size = 0;
	$main::last_verbose_state = 0;
      }

    #Return success
    return(0);
  }


##
## Subroutine that gets a line of input and accounts for carriage returns that
## many different platforms use instead of hard returns.  Note, it uses a
## global array reference variable ($infile_line_buffer) to keep track of
## buffered lines from multiple file handles.
##
sub getLine
  {
    my $file_handle = $_[0];

    #Set a global array variable if not already set
    $main::infile_line_buffer = {} if(!defined($main::infile_line_buffer));
    if(!exists($main::infile_line_buffer->{$file_handle}))
      {$main::infile_line_buffer->{$file_handle}->{FILE} = []}

    #If this sub was called in array context
    if(wantarray)
      {
	#Check to see if this file handle has anything remaining in its buffer
	#and if so return it with the rest
	if(scalar(@{$main::infile_line_buffer->{$file_handle}->{FILE}}) > 0)
	  {
	    return(@{$main::infile_line_buffer->{$file_handle}->{FILE}},
		   map
		   {
		     #If carriage returns were substituted and we haven't
		     #already issued a carriage return warning for this file
		     #handle
		     if(s/\r\n|\n\r|\r/\n/g &&
			!exists($main::infile_line_buffer->{$file_handle}
				->{WARNED}))
		       {
			 $main::infile_line_buffer->{$file_handle}->{WARNED}
			   = 1;
			 warning('Carriage returns were found in your file ',
				 'and replaced with hard returns.');
		       }
		     split(/(?<=\n)/,$_);
		   } <$file_handle>);
	  }
	
	#Otherwise return everything else
	return(map
	       {
		 #If carriage returns were substituted and we haven't already
		 #issued a carriage return warning for this file handle
		 if(s/\r\n|\n\r|\r/\n/g &&
		    !exists($main::infile_line_buffer->{$file_handle}
			    ->{WARNED}))
		   {
		     $main::infile_line_buffer->{$file_handle}->{WARNED}
		       = 1;
		     warning('Carriage returns were found in your file ',
			     'and replaced with hard returns.');
		   }
		 split(/(?<=\n)/,$_);
	       } <$file_handle>);
      }

    #If the file handle's buffer is empty, put more on
    if(scalar(@{$main::infile_line_buffer->{$file_handle}->{FILE}}) == 0)
      {
	my $line = <$file_handle>;
	#The following is to deal with files that have the eof character at the
	#end of the last line.  I may not have it completely right yet.
	if(defined($line))
	  {
	    if($line =~ s/\r\n|\n\r|\r/\n/g &&
	       !exists($main::infile_line_buffer->{$file_handle}->{WARNED}))
	      {
		$main::infile_line_buffer->{$file_handle}->{WARNED} = 1;
		warning('Carriage returns were found in your file and ',
			'replaced with hard returns.');
	      }
	    @{$main::infile_line_buffer->{$file_handle}->{FILE}} =
	      split(/(?<=\n)/,$line);
	  }
	else
	  {@{$main::infile_line_buffer->{$file_handle}->{FILE}} = ($line)}
      }

    #Shift off and return the first thing in the buffer for this file handle
    return($_ = shift(@{$main::infile_line_buffer->{$file_handle}->{FILE}}));
  }

##
## This subroutine allows the user to print debug messages containing the line
## of code where the debug print came from and a debug number.  Debug prints
## will only be printed (to STDERR) if the debug option is supplied on the
## command line.
##
sub debug
  {
    return(0) unless($DEBUG);

    $main::debug_number++;

    #Gather and concatenate the error message and split on hard returns
    my @debug_message = split("\n",join('',@_));
    pop(@debug_message) if($debug_message[-1] !~ /\S/);

    #Assign the values from the calling subroutine
    #but if called from main, assign the values from main
    my($junk1,$junk2,$line_num,$calling_sub);
    (($junk1,$junk2,$line_num,$calling_sub) = caller(1)) ||
      (($junk1,$junk2,$line_num) = caller());

    #Edit the calling subroutine string
    $calling_sub =~ s/^.*?::(.+)$/$1:/ if(defined($calling_sub));

    my $leader_string = "DEBUG$main::debug_number:LINE$line_num:" .
      (defined($calling_sub) ? $calling_sub : '') .
	' ';

    #Figure out the length of the first line of the error
    my $debug_length = length(($debug_message[0] =~ /\S/ ?
			       $leader_string : '') .
			      $debug_message[0]);

    #Put location information at the beginning of each line of the message
    foreach my $line (@debug_message)
      {print STDERR (($line =~ /\S/ ? $leader_string : ''),
		     $line,
		     ($verbose &&
		      defined($main::last_verbose_state) &&
		      $main::last_verbose_state ?
		      ' ' x ($main::last_verbose_size - $debug_length) : ''),
		     "\n")}

    #Reset the verbose states if verbose is true
    if($verbose)
      {
	$main::last_verbose_size = 0;
	$main::last_verbose_state = 0;
      }

    #Return success
    return(0);
  }


##
## This sub marks the time (which it pushes onto an array) and in scalar
## context returns the time since the last mark by default or supplied mark
## (optional) In array context, the time between all marks is always returned
## regardless of a supplied mark index
## A mark is not made if a mark index is supplied
## Uses a global time_marks array reference
##
sub markTime
  {
    #Record the time
    my $time = time();

    #Set a global array variable if not already set to contain (as the first
    #element) the time the program started (NOTE: "$^T" is a perl variable that
    #contains the start time of the script)
    $main::time_marks = [$^T] if(!defined($main::time_marks));

    #Read in the time mark index or set the default value
    my $mark_index = (defined($_[0]) ? $_[0] : -1);  #Optional Default: -1

    #Error check the time mark index sent in
    if($mark_index > (scalar(@$main::time_marks) - 1))
      {
	error("Supplied time mark index is larger than the size of the ",
	      "time_marks array.\nThe last mark will be set.");
	$mark_index = -1;
      }

    #Calculate the time since the time recorded at the time mark index
    my $time_since_mark = $time - $main::time_marks->[$mark_index];

    #Add the current time to the time marks array
    push(@$main::time_marks,$time)
      if(!defined($_[0]) || scalar(@$main::time_marks) == 0);

    #If called in array context, return time between all marks
    if(wantarray)
      {
	if(scalar(@$main::time_marks) > 1)
	  {return(map {$main::time_marks->[$_ - 1] - $main::time_marks->[$_]}
		  (1..(scalar(@$main::time_marks) - 1)))}
	else
	  {return(())}
      }

    #Return the time since the time recorded at the supplied time mark index
    return($time_since_mark);
  }

##
## This subroutine reconstructs the command entered on the command line
## (excluding standard input and output redirects).  The intended use for this
## subroutine is for when a user wants the output to contain the input command
## parameters in order to keep track of what parameters go with which output
## files.
##
sub getCommand
  {
    my $perl_path_flag = $_[0];
    my($command);

    #Determine the script name
    my $script = $0;
    $script =~ s/^.*\/([^\/]+)$/$1/;

    #Put quotes around any parameters containing un-escaped spaces or astericks
    my $arguments = [@$preserve_args];
    foreach my $arg (@$arguments)
      {if($arg =~ /(?<!\\)[\s\*]/ || $arg eq '')
	 {$arg = "'" . $arg . "'"}}

    #Determine the perl path used (dependent on the `which` unix built-in)
    if($perl_path_flag)
      {
	$command = `which $^X`;
	chomp($command);
	$command .= ' ';
      }

    #Build the original command
    $command .= join(' ',($0,@$arguments));

    #Note, this sub doesn't add any redirected files in or out

    return($command);
  }

##
## This subroutine checks to see if a parameter is a single file with spaces in
## the name before doing a glob (which would break up the single file name
## improperly).  The purpose is to allow the user to enter a single input file
## name using double quotes and un-escaped spaces as is expected to work with
## many programs which accept individual files as opposed to sets of files.  If
## the user wants to enter multiple files, it is assumed that space delimiting
## will prompt the user to realize they need to escape the spaces in the file
## names.
##
sub sglob
  {
    my $command_line_string = $_[0];
    return(-e $command_line_string ?
	   $command_line_string : glob($command_line_string));
  }


sub printVersion
  {
    my $script = $0;
    $script =~ s/^.*\/([^\/]+)$/$1/;
    print(($verbose ? "$script Version " : ''),
	  $software_version_number,
	  "\n");
    verbose("Generated using perl_script_template.pl\n",
	    "Version $template_version_number\n",
	    "Robert W. Leach\n",
	    "robleach\@lanl.gov\n",
	    "5/8/2006\n",
	    "Los Alamos National Laboratory\n",
	    "Copyright 2006");
    return(0);
  }

#This subroutine is a check to see if input is user-entered via a TTY (result is non-
#zero) or directed in (result is zero)
sub isStandardInputFromTerminal
  {return(-t STDIN || eof(STDIN))}

#This subroutine is a check to see if prints are going to a TTY.  Note, explicit prints
#to STDOUT when another output handle is selected are not considered and may defeat this
#subroutine.
sub isStandardOutputToTerminal
  {return(-t STDOUT && select() eq 'main::STDOUT')}
