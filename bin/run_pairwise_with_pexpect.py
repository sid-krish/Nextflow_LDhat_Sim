#!/usr/bin/env python
import sys

import pexpect

recom_tract_len = sys.argv[1]
sites_file = sys.argv[2]
locs_file = sys.argv[3]
lk_table = sys.argv[4]

child = pexpect.spawn(f"pairwise -seq {sites_file} -loc {locs_file} -lk {lk_table} -prefix pairwise_", encoding='utf-8',
                      logfile=sys.stdout)

child.expect_exact("Input average tract length for conversion model: ")
child.sendline(recom_tract_len)

child.expect_exact(
    "Do you wish to change grid over which to estimate likelihoods for (default = 101 points, 4Ner 0 - 100.0) (1/0) :")
child.sendline("0")

child.expect_exact("Do you wish to carry out a sliding windows analysis? (yes=1/no=0):")
child.sendline("0")

child.expect_exact("(0=No, 1=Total only, 2=Full table):")
child.sendline("0")

child.expect_exact("Estimate 4Ner by moment method? (yes=1, no=0)")
child.sendline("0")

child.expect_exact("Do you wish to test for recombination? (yes=1, no=0):")
child.sendline("0")

child.expect_exact(
    "Do you wish to test constant-rate model and estimate sampling distribution by simulation? (yes=1/no=0):")
child.sendline("0")

# from docs: You can also just expect the EOF if you are waiting for all output of a child to finish.
child.expect(pexpect.EOF)
