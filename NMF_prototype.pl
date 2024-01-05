# 
#   Prototype for Non-Negative Factorization 
#   Pablo Tamayo and Jean-Philippe Brunet,  May 6, 2003
#
#   The algorithm implemented is from:
#   D. Lee and H. S. Seung "Algorithms for Non-negative Factorization"
#   Advances in Neural Information Systems 13, 556 (2001)
#
#   See also Lee and Seung original article in Nature 401, 788 (1999).

BEGIN{
    open(STDERR, ">&STDOUT");
    $| = 1;                         # fflush() invoked at every output.
    print '';
}

use POSIX;
use Getopt::Std;
getopt('snfvcrpgmo');

$SIMULT_UPDATE = 0;    # do sequential update of H and W 
# $SIMULT_UPDATE = 1;    # do simultaneous update of H and W. DANGER: it violates assumption 
                        # that update of H is done with W fixed and viceversa and 
                        #  therefore the Euclidean error does not decrease monotonically 

$seed = 1234; # Default seed = 1234
$seed = $opt_s if defined($opt_s); # -s <random seed>

$niter = 500; # Default number of iterations = 500
$niter = $opt_n if defined($opt_n);  # -n <number of iterations>

$nfact = 4; # Default number of factors = 4
$nfact = $opt_f if defined($opt_f);  # -f <number of factors>

$v_file = "NULL";  # Default filename = NULL
$v_file = $opt_v if defined($opt_v);  # -v <input file>

$pflag = 0;  # Default print flag = 0 (false)
$pflag = $opt_p if defined($opt_p); # -p <print flag>

$norm_c = 0;  # Default normalize cols flag = 0 (false)
$norm_c = $opt_c if defined($opt_c); # -c <normalize columns flag>

$norm_r = 1; # Default normalize rows flag = 1 (true)
$norm_r = $opt_r if defined($opt_r); # -r <normalize rows flag>

$global_shift = 4;  # Default global shift = -1 (global shift such that minimum = 0
$global_shift = $opt_g if defined($opt_g); # -g <global shift>

$norm_m = 0; # Normalization type = 0 (default stardardize)
$norm_m = $opt_m if defined($opt_m); # -m <norm. type 0= stardardize, 1 = min=0, max=1>

$output_file = "NULL"; # Default output filename
$output_file = $opt_o if defined($opt_o); # -o <output file>

$filename = ">" . "$output_file";
open(OUT, $filename) || die "cannot open file: $filename \n";

print OUT "Running...\n";

srand($seed);

@v = ();
@vp = ();
@v_temp = ();
@vp_temp = ();
@w = ();
@h = ();
@h_old = ();
@v_sorted = ();
@vp_sorted = ();
@w_sorted = ();
@h_sorted = ();

@Euclid_error = ();

# Read input matrix V

open(V, "$v_file") || die "Cannot open V file: $v_file \n";

$nrows = 0;

print OUT "\n Reading input file (gct)...";

$colnames = ();
$rownames = ();
$rowdescs = ();

if ($v_file =~ /(.*).gct/) {
    $line = <V>;  # skip version line
    $line = <V>;  # skip header
    $line = <V>;  # read sample names
    chomp $line;
    @data = split("\t", $line);
    $datum = shift(@data); # skip accession
    $datum = shift(@data); # skip description
    $ncols = 0;
    while (defined($datum = shift(@data))) {
	$colnames[$ncols] = $datum;
	$ncols++;
    }

    while (<V>) {
	chomp $_;
	next if /^\s*$/;
	next if /^#/;
	@data = split("\t", $_);
	$rownames[$nrows] = shift(@data); # accession
	$rowdescs[$nrows] = shift(@data); # description
	$ncols = 0;
	while (defined($datum = shift(@data))) {
	    $v[$nrows][$ncols] = $datum;
#	print OUT "V elem, i=$nrows, j=$ncols, v=$v[$nrows][$ncols] \n";
	    $ncols++;
	}
	$nrows++;
	print OUT ".";
    }
} else {
    print OUT "\n Reading input file..."; 
    while (<V>) {
	chomp $_;
	next if /^\s*$/;
	next if /^#/;
#    @data = split("\t", $_);
	@data = split(" ", $_);
	$ncols = 0;
	while (defined($datum = shift(@data))) {
	    $v[$nrows][$ncols] = $datum;
#	print OUT "V elem, i=$nrows, j=$ncols, v=$v[$nrows][$ncols] \n";
	    $ncols++;
	}
	$nrows++;
	print OUT ".";
    }
}

print OUT "\n";

close V;

print OUT "\n total rows read: $nrows \n";
print OUT "\n total cols read: $ncols \n";

# Normalize columns if requested

if ($norm_c == 1) {
    print OUT "\n\n Normalizing columns...\n\n";
    for $j (0..$ncols - 1) { 
	$mean = 0;
	$sigma = 0;
	$max = $v[0][$j];
	$min = $v[0][$j];
	for $i (0..$nrows - 1) { 
	    $max = $v[$i][$j] if $v[$i][$j] > $max;
	    $min = $v[$i][$j] if $v[$i][$j] < $min;
	    $mean = $mean + $v[$i][$j];
	    $sigma = $sigma + $v[$i][$j]*$v[$i][$j];
	}
	$mean = $mean/$nrows;
	$sigma = $sigma/$nrows;
	$sigma = $sigma - $mean*$mean;
	$sigma = ($sigma > 0)? sqrt($sigma) : .001;
	for $i (0..$nrows - 1) { 
	    if ($norm_m == 0) {
		$v[$i][$j] = ($v[$i][$j] - $mean)/$sigma;
	    } else {
		if (($max - $min) > 0) {
		    $v[$i][$j] = ($v[$i][$j] - $min)/($max - $min);
		} else {
		    $v[$i][$j] = 0;
		}
	    }
	}
    }
}

# Normalize rows if requested

if ($norm_r == 1) {
    print OUT "\n\n Normalizing rows...\n\n";
    for $i (0..$nrows - 1) { 
	$mean = 0;
	$sigma = 0;
	$max = $v[$i][0];
	$min = $v[$i][0];
	for $j (0..$ncols - 1) { 
	    $max = $v[$i][$j] if $v[$i][$j] > $max;
	    $min = $v[$i][$j] if $v[$i][$j] < $min;
	    $mean = $mean + $v[$i][$j];
	    $sigma = $sigma + $v[$i][$j]*$v[$i][$j];
	}
	$mean = $mean/$ncols;
	$sigma = $sigma/$ncols;
	$sigma = $sigma - $mean*$mean;
	$sigma = ($sigma > 0)? sqrt($sigma) : .001;
	for $j (0..$ncols - 1) { 
	    if ($norm_m == 0) {
		$v[$i][$j] = ($v[$i][$j] - $mean)/$sigma;
	    } else {
		if (($max - $min) > 0) {
		    $v[$i][$j] = ($v[$i][$j] - $min)/($max - $min);
		} else {
		    $v[$i][$j] = 0;
		}
	    }
	}
    }
}

# Make global shift 

if ($global_shift > 0) {
    print OUT "\n\n Making global shift = $global_shift ...\n\n";
    for $i (0..$nrows - 1) { 
	for $j (0..$ncols - 1) { 
	    $v[$i][$j] = $global_shift + $v[$i][$j];
	}
    }
} elsif ($global_shift == -1) {
    print OUT "\n\n Making global shift so that minimum value is zero ...\n\n";
    $min = $v[0][0];
    for $i (0..$nrows - 1) { 
	for $j (0..$ncols - 1) { 
	    if ($v[$i][$j] < $min) {
		$min = $v[$i][$j];
	    }
	}
    }
    $global_shift = - $min;
    print OUT "\n\n Making global shift = $global_shift ...\n\n";
    for $i (0..$nrows - 1) { 
	for $j (0..$ncols - 1) { 
	    $v[$i][$j] = $global_shift + $v[$i][$j];
	}
    }
}

# Truncate negative values

print OUT "\n\n Truncating negative values...\n\n";
$warning = 0;
$nvals = 0;
for $i (0..$nrows - 1) {
   for $j (0..$ncols - 1) {
       if ($v[$i][$j] < 0.01) {
	   $v[$i][$j] = 0.01;
	   $nvals++;
	   if ($warning == 0) {
	       print OUT "\n WARNING: there is at least one zero or negative value at location ($i, $j) and potentially others \n\n";
	       $warning = 1;
	   }
       }
   }
}
print OUT "\n\n A total of $nvals negative values were truncated \n\n";

# initialize W and H with random numbers */

for $i (0..$nrows - 1) {
   for $j (0..$nfact - 1) {
       $w[$i][$j] = rand() + 0.01;
   }
}

for $i (0..$nfact - 1) {
    for $j (0..$ncols - 1) {
       $h[$i][$j] = rand() + 0.01;
   }
}


# NMF iteration loop

print OUT "\n Computing NMF \n";

  for $t (1..$niter) {
#    print OUT ".";

# Compute VP = W H   and Euclidean error =  ||V - W H||^2 

    $Euclid_error[$t] = 0.0;

    for $i (0..$nrows - 1) {
      for $j (0..$ncols - 1) {
	$vp[$i][$j] = 0.0;
        for $k (0..$nfact - 1) {
	  $vp[$i][$j] += $w[$i][$k] * $h[$k][$j];
	}
	$Euclid_error[$t] += ($v[$i][$j] - $vp[$i][$j])*($v[$i][$j] - $vp[$i][$j]);
      }
    }

    $Euclid_error[$t] = sqrt($Euclid_error[$t]/($nrows*$ncols));

    if ($t < 20) {
	$val = sprintf("%4.6f", $Euclid_error[$t]);
	print OUT "\n time = $t, Euclidean error = $val";
    } else {
      if ($t  % 20 == 0) {
	  $val = sprintf("%4.6f", $Euclid_error[$t]);
	  print OUT "\n time= $t,  Euclidean error = $val";
      }
    }

# Update H <- H (Wt V)/(Wt W H) */

    for $k  (0..$nfact - 1) {
      for $j (0..$ncols - 1) {
	$num_kj = 0.0;
	$den_kj = 0.0;
	for $i (0..$nrows - 1) {
	  $num_kj += $w[$i][$k] * $v[$i][$j];
	  $den_kj += $w[$i][$k] * $vp[$i][$j];
	}
	$h_old[$k][$j] = $h[$k][$j];
	$h[$k][$j] = $h[$k][$j]*$num_kj/$den_kj;
      }
    }

#  Update Vp 

    if ($SIMULT_UPDATE == 0) {
      for $i (0..$nrows - 1) {
	for $j (0..$ncols - 1) {
	  $vp[$i][$j] = 0.0;
	  for $k (0..$nfact - 1) {
	    $vp[$i][$j] += $w[$i][$k] * $h[$k][$j];
	  }
	}
      }
    }

# Update W <- W (V Ht)/(W H Ht) 
    
    for $i (0..$nrows - 1) {
      for $k (0..$nfact - 1) {
	$num_ik = 0.0;
	$den_ik = 0.0;
	for $j (0..$ncols - 1) {
	  if ($SIMULT_UPDATE == 0) {
	    $num_ik += $v[$i][$j] * $h[$k][$j];
	    $den_ik += $vp[$i][$j] * $h[$k][$j];
	  } else {
	    $num_ik += $v[$i][$j] * $h_old[$k][$j];
	    $den_ik += $vp[$i][$j] * $h_old[$k][$j];
	  }
	}
	$w[$i][$k] = $w[$i][$k]*$num_ik/$den_ik;
      }
    }

  }   # end of iteration loop

 print OUT "\n \n";


# cluster membership for H columns

@c_index = ();
@c_index_val = ();
for $i (0..$ncols - 1) {
    $c_index[$i] = 0;
    for $k (0..$nfact - 1) {
	if ($h[$k][$i] >= $h[$c_index[$i]][$i]) {
	    $c_index[$i] = $k;
	    $c_index_val[$i] = $h[$k][$i];
	}
    }
}

@c_count = ();
for $k (0..$nfact - 1) {
    $c_count[$k] = 0;
}
print OUT "\n\n";
for $i (0..$ncols - 1) {
    print OUT "\n case: $i factor: $c_index[$i]" if ($i < 100);
    $c_count[$c_index[$i]]++;
}
print OUT "\n\n";
$c_total = 0;
for $k (0..$nfact - 1) {
    print OUT "\n Factor: $k count: $c_count[$k]";
    $c_total = $c_total + $c_count[$k];
}

print OUT "\n\n (Col) total = $c_total";

# cluster membership for W rows

@r_index = ();
@r_index_val = ();
for $i (0..$nrows - 1) {
    $r_index[$i] = 0;
    for $k (0..$nfact - 1) {
	if ($w[$i][$k] >= $w[$i][$r_index[$i]]) {
	    $r_index[$i] = $k;
	    $r_index_val[$i] = $w[$i][$k];
	}
    }
}

@r_count = ();
for $k (0..$nfact - 1) {
    $r_count[$k] = 0;
}
print OUT "\n\n";
for $i (0..$nrows - 1) {
    print OUT "\n attribute: $i factor: $r_index[$i]" if ($i < 100);
    $r_count[$r_index[$i]]++;
}
print OUT "\n\n";
$r_total = 0;
for $k (0..$nfact - 1) {
    print OUT "\n Factor: $k count: $r_count[$k]";
    $r_total = $r_total + $r_count[$k];
}

print OUT "\n\n (Row) total = $r_total";

# Produce factor-ordered matrices

# Find column and row sort indices

%value = ();
@c_sort_order = ();
for $i (0..$ncols - 1) {
    $c_sort_order[$i] = $i;
    $value{$i} = $c_index[$i]*1000 + $c_index_val[$i];
}
@c_sort_order = sort byvalue @c_sort_order;
@colnames_s = ();
for $i (0..$ncols - 1) {
#    $colnames_s[$i] = "$colnames[$c_sort_order[$i]]" . "_" . "f$c_index[$c_sort_order[$i]]" . "_" . "$value{$c_sort_order[$i]}";
    $colnames_s[$i] = "$colnames[$c_sort_order[$i]]" . "_" . "f$c_index[$c_sort_order[$i]]";
}

%value = ();
@r_sort_order = ();
for $i (0..$nrows - 1) {
    $r_sort_order[$i] = $i;
    $value{$i} = $r_index[$i]*1000 + $r_index_val[$i];
}
@r_sort_order = sort byvalue @r_sort_order;
@rownames_s = ();
@rowdescs_s = ();
for $i (0..$nrows - 1) {
#    $rownames_s[$i] = "$rownames[$r_sort_order[$i]]" . "_" . "f$r_index[$r_sort_order[$i]]" . "_" . "$value{$r_sort_order[$i]}";
    $rownames_s[$i] = "$rownames[$r_sort_order[$i]]" . "_" . "f$r_index[$r_sort_order[$i]]";
    $rowdescs_s[$i] = "$rowdescs[$r_sort_order[$i]]";
}

# Sort H

for $i (0..$ncols - 1) {
    for $k (0..$nfact - 1) {
	$h_sorted[$k][$i] = $h[$k][$c_sort_order[$i]];
    }
}

# Sort W

for $i (0..$nrows - 1) {
    for $k (0..$nfact - 1) {
	$w_sorted[$i][$k] = $w[$r_sort_order[$i]][$k];
    }
}

# Sort columns and rows of V and VP (W H)

for $i (0..$nrows - 1) {
    for $j (0..$ncols - 1) {
	$v_sorted[$i][$j] = $v[$r_sort_order[$i]][$c_sort_order[$j]];
        $vp_sorted[$i][$j] = $vp[$r_sort_order[$i]][$c_sort_order[$j]];
    }
}

# Print (output) matrices as GCT files

print OUT "\n \n Done sorting matrices...\n";


print OUT "\n \n Printing matrices...\n\n";

$filename = ">" . "V.gct";
open(VGCT, $filename) || die "cannot open file: $filename \n";
$filename = ">" . "V_s.gct";
open(VGCTS, $filename) || die "cannot open file: $filename \n";
$filename = ">" . "W.gct";
open(WGCT, $filename) || die "cannot open file: $filename \n";
$filename = ">" . "H.gct";
open(HGCT, $filename) || die "cannot open file: $filename \n";
$filename = ">" . "WH.gct";
open(WHGCT, $filename) || die "cannot open file: $filename \n";
$filename = ">" . "WH_s.gct";
open(WHGCTS, $filename) || die "cannot open file: $filename \n";

$filename = ">" . "H_s.gct";
open(HGCTS, $filename) || die "cannot open file: $filename \n";
$filename = ">" . "W_s.gct";
open(WGCTS, $filename) || die "cannot open file: $filename \n";

print VGCT "#1.2\n";
print VGCTS "#1.2\n";
print WGCT "#1.2\n";
print WGCTS "#1.2\n";
print HGCT "#1.2\n";
print HGCTS "#1.2\n";
print WHGCT "#1.2\n";
print WHGCTS "#1.2\n";

print VGCT "$nrows\t$ncols\n";    
print VGCTS "$nrows\t$ncols\n";
print WGCT "$nrows\t$nfact\n";
print WGCTS "$nrows\t$nfact\n";
print HGCT "$nfact\t$ncols\n";
print HGCTS "$nfact\t$ncols\n";
print WHGCT "$nrows\t$ncols\n";
print WHGCTS "$nrows\t$ncols\n";

print VGCT "Name\tDescription\t";
print VGCTS "Name\tDescription\t";
print WGCT "Name\tDescription\t";    
print WGCTS "Name\tDescription\t";    
print HGCT "Name\tDescription\t";
print HGCTS "Name\tDescription\t";
print WHGCTS "Name\tDescription\t";

for $i (0..$ncols - 1) {
    print VGCT "$colnames[$i]\t";
    print VGCTS "$colnames_s[$i]\t";
    print HGCT "$colnames[$i]\t";
    print HGCTS "$colnames_s[$i]\t";
    print WHGCT "$colnames[$i]\t";
    print WHGCTS "$colnames_s[$i]\t";
}
print VGCT "\n";
print VGCTS "\n";
print HGCT "\n";
print HGCTS "\n";
print WHGCT "\n";
print WHGCTS "\n";

for $i (0..$nfact - 1) {
    print WGCT "f$i\t";
    print WGCTS "f$i\t";
}
print WGCT "\n";
print WGCTS "\n";

print OUT "\n Input V matrix from file $v_file after normalization, global shift and truncation  \n" if $pflag == 1;
for $i (0..$nrows - 1) {
    print VGCT "$rownames[$i]\t$rowdescs[$i]\t";
    for $j (0..$ncols - 1) {
	$val = sprintf("%5.3f", $v[$i][$j]);
	print OUT "$val " if $pflag == 1;
	print VGCT "$val\t";
    }
    print OUT "\n" if $pflag == 1;
    print VGCT "\n";
}
print OUT "\n " if $pflag == 1;

print OUT "\n Sorted V  matrix \n" if $pflag == 1;
for $i (0..$nrows - 1) {
    print VGCTS "$rownames_s[$i]\t$rowdescs_s[$i]\t";
    for $j (0..$ncols - 1) {
	$val = sprintf("%5.3f", $v_sorted[$i][$j]);
	print OUT "$val " if $pflag == 1;
	print VGCTS "$val\t";
    }
    print OUT "\n" if $pflag == 1;
    print VGCTS "\n";
}
print OUT "\n " if $pflag == 1;

print OUT "\n W matrix \n" if $pflag == 1;
for $i (0..$nrows - 1) {
    print WGCT "$rownames[$i]\t$rowdescs[$i]\t";
    for $k (0..$nfact - 1) {
	$val = sprintf("%5.3f", $w[$i][$k]);
	print OUT "$val " if $pflag == 1;
	print WGCT "$val\t";
    }
    print OUT "\n" if $pflag == 1;
    print WGCT "\n";

}
print OUT "\n" if $pflag == 1;

print OUT "\n Sorted W matrix \n" if $pflag == 1;
for $i (0..$nrows - 1) {
    print WGCTS "$rownames_s[$i]\t$rowdescs_s[$i]\t";
    for $k (0..$nfact - 1) {
	$val = sprintf("%5.3f", $w_sorted[$i][$k]);
	print OUT "$val " if $pflag == 1;
	print WGCTS "$val\t";
    }
    print OUT "\n" if $pflag == 1;
    print WGCTS "\n";

}
print OUT "\n" if $pflag == 1;

print OUT "\n H matrix \n" if $pflag == 1;
for $k (0..$nfact - 1) {
    print HGCT "f$k\tf$k\t";
    for $j (0..$ncols - 1) {
	$val = sprintf("%5.3f", $h[$k][$j]);
	print OUT "$val " if $pflag == 1;
	print HGCT "$val\t";
    }
    print OUT "\n" if $pflag == 1;
    print HGCT "\n";
}
print OUT "\n " if $pflag == 1;

print OUT "\n Sorted H matrix \n" if $pflag == 1;
for $k (0..$nfact - 1) {
    print HGCTS "f$k\tf$k\t";
    for $j (0..$ncols - 1) {
	$val = sprintf("%5.3f", $h_sorted[$k][$j]);
	print OUT "$val " if $pflag == 1;
	print HGCTS "$val\t";
    }
    print OUT "\n" if $pflag == 1;
    print HGCTS "\n";
}
print OUT "\n " if $pflag == 1;

print OUT "\n Vp = W H  matrix \n" if $pflag == 1;
for $i (0..$nrows - 1) {
    print WHGCT "$rownames[$i]\t$rowdescs[$i]\t";
    for $j (0..$ncols - 1) {
	$val = sprintf("%5.3f", $vp[$i][$j]);
	print OUT "$val " if $pflag == 1;
	print WHGCT "$val\t";
    }
    print OUT "\n" if $pflag == 1;
    print WHGCT "\n";
}
print OUT "\n " if $pflag == 1;

print OUT "\n Sorted Vp = W H  matrix \n" if $pflag == 1;
for $i (0..$nrows - 1) {
    print WHGCTS "$rownames_s[$i]\t$rowdescs_s[$i]\t";
    for $j (0..$ncols - 1) {
	$val = sprintf("%5.3f", $vp_sorted[$i][$j]);
	print OUT "$val " if $pflag == 1;
	print WHGCTS "$val\t";
    }
    print OUT "\n" if $pflag == 1;
    print WHGCTS "\n";
}
print OUT "\n " if $pflag == 1;

print OUT "\n Error matrix  V - W H \n" if $pflag == 1;
for $i (0..$nrows - 1) {
    for $j (0..$ncols - 1) {
	$val = sprintf("%5.3f", $v[$i][$j] - $vp[$i][$j]);
	print OUT "$val " if $pflag == 1;
    }
    print OUT "\n" if $pflag == 1;
}
print OUT "\n " if $pflag == 1;

close OUT;

sub byvalue {
    $value{$a} <=> $value{$b};
}
