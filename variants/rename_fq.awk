#!/bin/awk -f

FNR == 1 {
    print "Analyzing file", FILENAME > "/dev/tty"
    n = split(FILENAME, filename_split, "/")
    basename = filename_split[n]
    m = split(basename, basename_split, ".")
    suffix = basename_split[1]"."basename_split[2]
}
{
    if (NR % 4 == 1)
        print $0":"suffix
    else
        print $0
}
