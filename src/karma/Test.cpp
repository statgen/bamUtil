/*
 * Copyright (c) 2009 Regents of the University of Michigan
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "MapperSE.h"
#include "MapperSEBaseSpace.h"
#include "MapperSEColorSpace.h"
#include "MapperPE.h"
#include "MapperPEBaseSpace.h"
#include "Test.h"

#include <iomanip>
#include <iostream>
#include <fstream>

int Test::testSEBaseSpaceReads()
{
    MapperUserOptions mapperOptions;

    MapperSEBaseSpace se;
    std::string readFrag;
    std::string qualityFrag;
    std::string tag = "";
    int totalTests = 0;
    int success = 0;

    mapperOptions.mismatchCutoff = 4;
    mapperOptions.debug = false;

    se.initMapper(gs, wordIndex, wordHashLeft, wordHashRight, mapperOptions);

#if 0
    //
    // This read, although it is 35 bases, has three
    // exact matches in the genome.  It was taken from
    // chromosome 19, index 46374413, but maps to two
    // other locations.
    //
    // I'll leave it here, but commented out until we
    // improve the test to handle numBestMatch.
    //
    // @Chromosome_19_046374413_Genome_2691223814
    // GACAGAGTCTTGCTCTCTCACCCAGGCTGAAGTGC
    // +
    // 55555555555555555555555555555555555
    //
    success += !se.test(++totalTests,                   // test number
                        "GACAGAGTCTTGCTCTCTCACCCAGGCTGAAGTGC",   // read
                        "55555555555555555555555555555555555",   // quality
                        'F',                                     // direction
                        19,                                      // chr
                        46374413,                                // genomeIndex_t
                        0,                                       // misMatches
                        0);                                      // quality
#endif

    //
    // from /home/zhanxw/compareMapSoft/SEdbSNP/single1.fastq
    // @Chromosome_21_010047358_Genome_2781144374
    // CCAGTCGAAGACATTACAAGATCGTAAGACTACAG READ
    // CCAGACGAAGACATTACAAGAATGTAAGACTACAG REFERENCE
    //     ^                ^^
    // This is not a terribly special read... just three mismatches...
    //
    success += !se.test(++totalTests,                   // test number
                        "CCAGTCGAAGACATTACAAGATCGTAAGACTACAG",   // read
                        "55555555555555555555555555555555555",   // quality
                        'F',                                     // direction
                        21,                                      // chr
                        10047358,                                // genomeIndex_t
                        3,                                       // misMatches
                        60,                                      // quality
                        "35M");                                  // CIGAR

    // first word appears 3X in whole genome, 2nd word is high repeat
    // (=0x2fffffff)
    // For faster testing with this read:
    // head -c 33000000 NCBI36.fa >FOO.fa
    // ./karma --createIndex --reference FOO.fa --maxRepeat 400
    // Chromosome_1__027449926_Genome_0027449925       4       1       27449926        0       35M             0       0       GCTGTTAATTGACTTGTTTTTTTTTTTTTTCTTTT     55555555555555555555555555555555555     HA:i:3  ER:Z:repeats
    //
    success += !se.test(++totalTests,                   // test number
                        "GCTGTTAATTGACTTGTTTTTTTTTTTTTTCTTTT",   // read
                        "55555555555555555555555555555555555",   // quality
                        'F',                                     // direction
                        1,                                       // chr
                        27449926,                                // genomeIndex_t
                        0,                                       // misMatches
                        0,                                       // quality
                        "35M");                                  // CIGAR

    // this test is handy because there is only a single word that matches:
    // NOTE:
    // the real sequence is 1 bp difference than the following test read
    // reference seq: CCCGG AACGA TACCA AAGTC GCGGG CTTCT TCCGCCT
    //   test string: CCCGG AACGA TACCA AAGTC GCGGG CTTCC TCCGCCT
    //                                                  ^ different.
    // However, this mismatch will not count,the reason is by default,
    // readIndexCutoff = 24;   // maq uses this by default
    // and mismatch only counts the number of mismatches within the first
    // readIndexCutoff bps.
    success += !se.test(++totalTests,                   // test number
                        "CCCGGAACGATACCAAAGTCGCGGGCTTCCTCCGCCT", // read
                        "7AB8A<'BB+9;B;36;17A<B<B&8(>4'0<95:2&", // quality
                        'F',                                     // direction
                        5,                                       // chr
                        176762769,                               // genomeIndex_t
                        1,                                       // misMatches
                        6,                                       // quality
                        "37M");                                  // CIGAR
    success += !se.test(++totalTests,                   // test number
                        "GCAAAAGAACCCAATTTAACTTTACTAAACCAAATG",  // read
                        "@@@AAAAB9BB2;6BBBBB:BA8BBBBB@@;6@@(7",  // quality
                        'F',                     // direction
                        1,                         // chr
                        34994690,                     // chromosomeIndex
                        2,                         // misMatches
                        47,                         // quality
                        "36M");                                  // CIGAR

    //
    // from Xiaowei: HWI-EAS159:8:108:1142:1296#0/1
    //
    // This is test helped find a bug where the match position was
    // smaller than 15, so we underflowed doing a subtract in
    // MapperBase::evalAllCandidatePositions
    // The word starting at read index 0x17057057 maps to the
    // genome at index 4, 10, 16, 22, ... so is a classis
    // boundary condition problem.
    //
    //
    success += !se.test(++totalTests,                   // test number
                        "AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTACCCCTAACCC",  // read
                        "b_Ya`Zbb_aa_aa^aa`a___a]]_`a__a___Y^`__\\ZYX^\\[ZV[RBBBBBBBBBBBBBBBBBBBBB",  // quality
                        'F',                     // direction
                        5,                         // chr
                        63908,                     // chromosomeIndex
                        1,                         // misMatches
                        33,                         // quality
                        "71M");                  // CIGAR

#if 0
    // possibly useful for later.
    // this appears to me to be a candidate for SW alignment - pha
    //
    // from /home0/ftp.1000genomes.ebi.ac.uk/data/NA18000/sequence_read/SRR013948_1.recal.fastq.gz, tag @SRR013948.484 30PT5AAXX:5:1:4:372 length=76:
    //
    success += !se.test(++totalTests,
                        "CTCCAGCCTGGATGACAGAGCTAGACTCANCCCCCNNNNNNNNNNNNNNNCCCCCGCCNCCCCNCNNCCCCNCNNN",
                        "??>7@?>@@8@6D:585?<C4CC?63D)(%AAD=F%&&#('%%&%'%'#&@A?=B%>@%A<:C&>&&B8:D%=%&&",
                        'F',
                        1,
                        164084460,
                        1,
                        428);
#endif

    // this test use the reverse complement read of test 1
    success += !se.test(++totalTests,                   // test number
                        "AGGCGGAGGAAGCCCGCGACTTTGGTATCGTTCCGGG", // read
                        "&2:59<0'4>(8&B<B<A71;63;B;9+BB'<A8BA7", // quality
                        'R',                                     // direction
                        5,                                       // chr
                        176762769,                               // chromosomeIndex
                        1,                                       // misMatches
                        6,                                       // quality
                        "37M");                                  // CIGAR


    // @Chromosome_5__146523451_Genome_1027499208
    //
    // This is a special case - all the index words point to high repeat
    // index words, since this is now handled using the secondary hash,
    // it should map just fine.
    //
    success += !se.test(++totalTests,                   // test number
                        "CAGCATTTTGGGAGGCCGAGGCAGGTGGATCATGAGGTCAGGAGATCAAGACCATCCTGGCTAACACAGTGAAACCCTGTCTCTACTAAAAACACAAAAA",
                        "5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555",
                        'F',
                        5,
                        146523451,
                        0,
                        0,
                        "100M");                                   // CIGAR

    //
    // sanger tests from a broken PE test:
    //
    // This is a synthetic read known to originate from the
    // given location.
    //
    // This read is tricky because it has 12 mismatches,
    // but a reasonably decent sumQ of 66.
    //
    // A simple arbitrary mismatch count cutoff will eliminate this
    // read from that position.
    //
    success += !se.test(++totalTests,                   // test number
                        "AACCAAAATCATTTGAATAGGCTGAGAGTGAGACAAATAAAACAAATGCCCCACAATGAGAGAATGATACTGAAAACCCCCAAGGATGAATGCGATCTCATCCATCAC",
                        "%4======<<8<<=;63<8:95:=95.:6//'--.%%80)0))0;2/;4)%%:3%%+%20)%;2)+4<2%%%+02%)8.%2-22+4))17+)%(/-%1%%%1..%%%1",
                        'F',
                        9,                                          // chromosome
                        83445275,                                   // chromosomeIndex
                        12,                                         // misMatches
                        66,                                         // quality
                        "108M");                                    // CIGAR

    // this read is the pair of the above read and cleanly maps:
    success += !se.test(++totalTests,                   // test number
                        "ACTATATCATTACAGTAACACTCTGAATGTAAACAGCAAAGTCAGTTTCATTATAGTAATATAATTGACCCAAGCTGCACTAAGTAAAATGTGAATTAAAATTCTCAA",
                        "%1==========<=7:;<:<=;==<<<==8,7<5-::;;85;9;;;9<<9:;9;9<89::/<=7(,5:96//13699661+1/953653%-/17-)))/24,7,%-'0",
                        'R',
                        9,                                          // chromosome
                        83445275+124,                               // chromosomeIndex
                        1,                                          // misMatches
                        28,                                         // quality
                        "108M");                                    // CIGAR

    mapperOptions.forceSmithWaterman = true;
    se.initMapper(gs, wordIndex, wordHashLeft, wordHashRight, mapperOptions);
    //
    // This is not an especially good test - we don't know anything about the validity of
    // the result, nor do we know the CIGAR string, although it should be something like #M.
    //
    success += !se.test(++totalTests,                   // test number
                        "TTAAAATCCTTCACCGAGCTTTTCACTTGGGAAAGG",
                        ">>?>>>???@@@@@@@/@<@AA;=2<<=/1.117;,",
                        'F',
                        1,                                          // chromosome
                        102675073,                                  // chromosomeIndex
                        0,                                          // misMatches
                        0,                                          // quality
                        "36M");                                     // CIGAR

    //
    // @Chromosome_8__098604584_Genome_1490159623
    // CTCCACCTCCCAGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGTGCCCACCATGGC
    //
    // no indels, but SW is forced on, so we expect to see same location
    //
    success += !se.test(++totalTests,                   // test number
                        "CTCCACCTCCCAGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGTGCCCACCATGGC",
                        "555555555555555555555555555555555555555555555555555555555555555555555555555",
                        'F',
                        8,                                          // chromosome
                        98604584,                                   // chromosomeIndex
                        0,                                          // misMatches
                        0,                                          // quality
                        "75M");                                     // CIGAR
#if 0
    //
    // @Chromosome_8__098604584_Genome_1490159623
    // CTCCACCTCCCAGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGTGCCCACCATGGC
    //            ^ delete this base
    //
    // currently, the chr index is wrong... off by 1
    //
    success += !se.test(++totalTests,                   // test number
                        "CTCCACCTCCCGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGTGCCCACCATGGC",
                        "55555555555555555555555555555555555555555555555555555555555555555555555555",
                        'F',
                        8,                                          // chromosome
                        98604584,                                   // chromosomeIndex
                        0,                                          // misMatches
                        50,                                         // quality
                        "11M1D63M");                                // CIGAR

    //
    // from /home4/ftp.1000genomes.ebi.ac.uk/ftp/pilot_data/technical/working/20090924_pilot3_bc_bams/NA12878.454.BCM.MOSAIK.SRP000033.2009_09.bam
    // read NA12878-SRR005999.528968
    //
    // MOSAIK maps this as:
    //
    // chr 1       index 695700  qual 46      CIGAR 82M1D345M1D9M2D46M
    //
    // XXX NOTE: KARMA gets a CIGAR string of "81M1D343M1D6M2D52M"
    // This is the same CIGAR as the MOSAIK one.  The ambiguity is
    // because of identical bases as the deletion positions, which
    // means that either one could be called as "deleted".
    //
    // bandsize has to be at least 5.
    //
    mapperOptions.smithWatermanBandSize = 8;
    se.setMapperOptions(mapperOptions);
    success += !se.test(++totalTests,                   // test number
                        "GCCTAGCTAATCTTTGTATTTCTAGTAGAGATGCGGTTTTGCCACATTGCCAGTCCACCTGTCTCAGCCCCGCAAAGTGCTGTATTACAGGAGTGAGCCACTGCACCCAGCATTTGCCAAGATCTTTGATGGCAGGCTTTTTCCAGGTGATCAGTCCTTGTCTGGTCTGGCTCTGCCCCACTCTCCTTCTCACCTAGTTGGAATCCCTAGCTACTTTTCAGTAGAGGAGAGTGTGTACCCCAATCCCAGCTTGGTTCAGATCTGCATTTAACTCATGGAACCTGGCTGCTCCCCAGGTCCTGAAGAAAAAAACGGTCTCTCTGTGGGTATGATAAAGGATGGGCCTGTCCCCAGGACCCTGTGAGAGGGAAGCCCAATGTCCCACCAGGTTGGCAGGGCTGGGGAAGGGAAAGTGTTATGGCAGCCCAAGAAAAAAGAGGCAGCAGAGGGAGCAGGACAGCGCTCACATGGAACTCATGCCA",
                        "::::::::::;::;;;;;;;;;;;;::::::::::::::::::::::::::::::::::::::::::;;;;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;;;;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;;;;;;:;;;;;;;;;;::;;;;;;;::::::;;:;;;;;:::::;;;;;;::::;:;:::::::;;::;;;;::::::::;;::::::::;;;:;;;;;;;;;;:;;;;:;;;;;;;;;;;;;;:;;;::::;;;;;;;;9999:9;;:::::::::::::999::::;;;::::;::;;:;;;;;;;::::9999997779999;;;;;:;:;::999;;;9999999599;;;;::::::;:;;::;;;;;;;;::;;;;::::::;;;;;:",
                        'F',
                        1,                                          // chromosome
                        695700,                                     // chromosomeIndex
                        0,                                          // misMatches
                        250,                                        // quality
                        "81M1D343M1D6M2D52M");                      // CIGAR
#endif
    std::cout << "Single end tests: " << success << "/" << totalTests << " pass."<< std::endl;
    return 0;
}

int Test::testSEColorSpaceReads()
{
    MapperUserOptions mapperOptions;

    MapperSEColorSpace se;
    std::string readFrag;
    std::string qualityFrag;
    std::string tag = "";
    int totalTests = 0;
    int success = totalTests;

    mapperOptions.mismatchCutoff = 4;
    mapperOptions.debug = false;

    se.initMapper(gs, wordIndex, wordHashLeft, wordHashRight, mapperOptions);

    // this test is handy because there is only a single word that matches:
    // and this reads is converted from single end base space read in test 1
    // we expect the match position to be 176762769+1
    // (since color space reference genome is starting from the second base pair of chr 1.)
    // also, the quality is 30 ( but I don't quite sure how to get 30)
    // note:
    // the reference is: 003020132331010021233300320220203302
    // the read is :     003020132331010021233300320202203302
    success += se.test(++totalTests,                                       // test number
                       "003020132331010021233300320202203302", // read
                       "7AB8A<'BB+9;B;36;17A<B<B&8(>4'0<95:2", // quality
                       'F',                                     // direction
                       5,                                       // chr
                       176762770,                               // chromosomeIndex
                       2,                                       // misMatches
                       30);                                      // quality
    // this test is the reverse complement of the last test
    success += se.test(++totalTests,                                       // test number
                       "203302202023003332120010133231020300", // read
                       "2:59<0'4>(8&B<B<A71;63;B;9+BB'<A8BA7", // quality
                       'R',                                     // direction
                       5,                                       // chr
                       176762770,                               // chromosomeIndex
                       2,                                       // misMatches
                       30);                                      // quality
    //
    // from: /home0/ftp.1000genomes.ebi.ac.uk/data/NA19240/sequence_read/5/SRR007215.fastq.gz
    // read tag: VAB_SOLiD0097_20080602_1_Pilot2_YRI_D_library_MP_1_1_8Kb1_36_1338
    //
    success += se.test(++totalTests,                                       // test number
                       "G3320030011111131213300222", // read
                       "!')'.04,78<7/48=&)$)-49%1*", // quality
                       'R',                                     // direction
                       5,                                       // chr
                       176762770,                               // chromosomeIndex
                       2,                                       // misMatches
                       30);                                      // quality
    std::cout << "Single end color space tests: " << success << "/" << totalTests << " pass."<< std::endl;
    return 0;
}


int Test::testPEBaseSpaceReads()
{
    MapperUserOptions mapperOptions;

    MapperPEBaseSpace pe1, pe2;
    std::string readFrag;
    std::string qualityFrag;
    std::string tag = "";
    int success = 0;
    int totalTests = 0;

    mapperOptions.mismatchCutoff = -1;
    mapperOptions.debug = false;

    pe1.initMapper(gs, wordIndex, wordHashLeft, wordHashRight, mapperOptions);
    pe2.initMapper(gs, wordIndex, wordHashLeft, wordHashRight, mapperOptions);

    //
    // NB: some of the quality strings have trigraph character sequences in them,
    // for example: ??< is equivalent to { --- this forces us to require -Wno-trigraphs with g++. 
    //              ??> is equivalent to } --- this forces us to require -Wno-trigraphs with g++. 
    success += !pe1.test(++totalTests,
                         &pe2,
                         "GAAAGGGGGCTATTCCTAGTTTTATTGCTATAGCCA", "CCCCDCDCEDDDC>CDDEDECBEBEDEC@ABD=BAA", 'R', 1, 555181, 0,
                         "CAACCGCATCCATAATCCTTCTAATAGCTATCCTCT", "@@?DAA@AADFFBGDADADAECBCAFBA?A@?=>??", 'F', 1, 555073, 0
                        );

    //
    // I used extractGenomeSequnece to confirm that both of the following reads are on
    // the forward strand.
    // numBest == 2, so 2nd read is random...
    success += !pe1.test(++totalTests,
                         &pe2,
                         "ATTCCATTCCATTCCCcTGtACTCGGGTTGATTCCA", ">>?>>>???@@@@@@@@@@@>A?@?9=9=90>>>:,", 'F', 10, 41699896, 2,
                         "TTgCATTCCATTCCATTCCATTCCATTCCATTCCAT", "354?87:<>;79?==99;?A;;>@;87<=<76???<", 'F', 10, 41700502, 1
                        );

    // note: the following pair has 2 best matches, so depending on
    // which the random number generator says to write, we'll get
    // one or the other (therefore the prior tests matter).
    //
    // I used extractGenomeSequnece to confirm that both of the following reads are on
    // the forward strand.
    //
    globalRandom.Reset(2);  // impirical ... different seeds may fail
    success += !pe1.test(++totalTests,
                         &pe2,
                         "CcTTaCTcTCCATTaCATTCCATTCCATTCGGGTTG", "@=??9???@@?>@@=A7@A@@@A?>>3?;C<:;>;;", 'F', 4, 48796294, 4,
                         "ATTCCATTCCATTCCATTCCATTCCATTACAaTCgA", "=4:>>=4:@AA<@?AA@@@@@@@@@@@???>>>?>>", 'F', 4, 48795762, 2
                        );


    //
    // from sanger synthesized dataset 12/22/2009 108mer_hap1_1_20091222.fastq.gz
    // 9:83445275:F:124;9,83445326,A,G;None/1 and
    // 108mer_hap1_2_20091222.fastq.gz:
    // 9:83445275:F:124;9,83445326,A,G;None/2
    //
    success += !pe1.test(++totalTests,
                         &pe2,
                         "AACCAAAATCATTTGAATAGGCTGAGAGTGAGACAAATAAAACAAATGCCCCACAATGAGAGAATGATACTGAAAACCCCCAAGGATGAATGCGATCTCATCCATCAC", "%4======<<8<<=;63<8:95:=95.:6//'--.%%80)0))0;2/;4)%%:3%%+%20)%;2)+4<2%%%+02%)8.%2-22+4))17+)%(/-%1%%%1..%%%1", 'F', 9, 83445275, 2,
                         "ACTATATCATTACAGTAACACTCTGAATGTAAACAGCAAAGTCAGTTTCATTATAGTAATATAATTGACCCAAGCTGCACTAAGTAAAATGTGAATTAAAATTCTCAA", "%1==========<=7:;<:<=;==<<<==8,7<5-::;;85;9;;;9<<9:;9;9<89::/<=7(,5:96//13699661+1/953653%-/17-)))/24,7,%-'0", 'R', 9, 83445275-124, 1
                        );


    // also from the above sanger dataset
    success += !pe1.test(++totalTests,
                         &pe2,
                         "GCACCATCCTTTCTAGCTGTACAGCCCTCCTGAGCCTTAGTTCTTACATCTTGAAGATGGAACCAGCTCAACAAGCATAGGGATGTAGCAAGAATCAAGACAATGTAG", "%0<0<8:/::<;<<:<<=<:9:<:4-4:<:;2719686:8:9236;90'.;0.1.0.3(7-5'+)2-86/8888;/:2)3<:66<9/(--'.)/)%%%%%))%)/+2%", 'F', 23, 55351604, 2,
                         "GTCAAGGTCATAGAAGAAATAGTGCCACACAAGGTAGAGGGTGAAGATAACTCTAGAGAACCATTCTCAAACCAACGGCTAGGCGCTGCATCCCGAGCGTATTGATTC", "%/>>>8<:><>>;<=<;2;<==>>;>><==<<=>;;==</+%%5,0313%/%%&%(+%%1+/8%%)%))/%%,+%).))&%+%+%&+%,,/.%)&0,+-)),%)-%-%", 'R', 23, 55351604-91, 1
                        );

    //
    // small insert from  alignment_comp_250110/108mer_hap2_[12].fastq.gz:
    //
    // @X:46303310:F:93;X,46303397,*,+6TTGTAG;None/1
    // CTCTTAAAAAAAAAATACTAAAACTAAGTTCTATTACAATTTATGTTTTACTTCTCATAAAATAAATTATACATATTATCAGAAGGAATTGTAGAAAAAATGTGTAAA
    // +
    // CC=CCCCC=DCBBCBACCB?<DBCCCAC;BBCCAB?BB=@@@=>:@@A=@>>BA9?=;7<C@@69@8?>@=:=?@B=C=60(4:5964=9?:8027;>8)9?1/84-'
    // and
    // @X:46303310:F:93;X,46303397,*,+6TTGTAG;None/2
    // ACTTACCAGGTCCTTTTGAGTTTTCTAATAAAAACATAAATCAATAAACAAACATATTTTCTGGGCAAAGAGTGGCTTCTCACAGGAATTTTTTCTAATTTTACACAT
    // +
    // BBBBBB@BBB>BBBBBBB>AB<40BBBBB@@@ABBABBBB?BBBBA@@>BBA723;@@=7:4+>@@@?;==8@68,8-+6(<@>8624,=B@B7'8@=(,7,0>3->B
    //
    // XXX this aligns now because the insert is near the end of the read
    //
    success += !pe1.test(++totalTests,
                         &pe2,
                         "CTCTTAAAAAAAAAATACTAAAACTAAGTTCTATTACAATTTATGTTTTACTTCTCATAAAATAAATTATACATATTATCAGAAGGAATTGTAGAAAAAATGTGTAAA", "CC=CCCCC=DCBBCBACCB?<DBCCCAC;BBCCAB?BB=@@@=>:@@A=@>>BA9?=;7<C@@69@8?>@=:=?@B=C=60(4:5964=9?:8027;>8)9?1/84-'", 'F', 23, 46303310, 16,
                         "ACTTACCAGGTCCTTTTGAGTTTTCTAATAAAAACATAAATCAATAAACAAACATATTTTCTGGGCAAAGAGTGGCTTCTCACAGGAATTTTTTCTAATTTTACACAT", "BBBBBB@BBB>BBBBBBB>AB<40BBBBB@@@ABBABBBB?BBBBA@@>BBA723;@@=7:4+>@@@?;==8@68,8-+6(<@>8624,=B@B7'8@=(,7,0>3->B", 'R', 23, 46303403, 0
                        );


    std::cout << "Paired end tests: " << success << "/" << totalTests << " pass."<< std::endl;
    return 0;
}


int Test::testRemapReference(std::string &outputFile, std::string &whichChromosome, int numberBases, int skipBases)
{
    std::ostream    *out = &std::cout;
    std::ofstream   outFile;

    MapperSEBaseSpace se;
    MapperUserOptions mapperOptions;

    mapperOptions.mismatchCutoff = 4;
    mapperOptions.debug = false;

    se.initMapper(gs, wordIndex, wordHashLeft, wordHashRight, mapperOptions);

    std::string quality = "";

    for (int j=0; j<numberBases; j++) quality += (char)('!' + 20);        // phred quality score 20

    std::string read;

    int firstChromosome = 0;
    int lastChromosome = gs->getChromosomeCount() - 1;

    if (outputFile != "")
    {
        outFile.open(outputFile.c_str(), std::ios_base::out | std::ios_base::trunc);
        out = &outFile;
    }
    if (whichChromosome != "")
    {
        int which = gs->getChromosome(whichChromosome.c_str());
        if (which < 0 || which >= gs->getChromosomeCount())
        {
            std::cerr << "The given reference does not have a chromosome named " << whichChromosome << std::endl;
            return 1;
        }
        firstChromosome = lastChromosome = which;
    }

    for (int chromosome = firstChromosome; chromosome <= lastChromosome; chromosome++)
    {

        genomeIndex_t   chromosomeStart = gs->getChromosomeStart(chromosome);

        for (genomeIndex_t chromosomePosition = 0;
                chromosomePosition < gs->getChromosomeSize(chromosome) - numberBases;
                chromosomePosition += skipBases)
        {

            genomeIndex_t readGenomePosition = chromosomeStart + chromosomePosition;

            gs->getString(read, chromosome, chromosomePosition + 1, numberBases);

            (*out) << gs->getChromosomeName(chromosome) << ":" << (chromosomePosition+1) << "\t";

            if (se.setReadAndQuality(read.c_str(), read.size(), quality.c_str()))
            {
                (*out) << "N" << std::endl;     // N -> too many N's or other bad bases to get 2 index words
                continue;
            }
            se.fragmentTag = "";
            se.MapSingleRead();

            if (se.getBestMatch().qualityIsValid())
            {
                int qualityScore = se.getBestMatch().getQualityScore();

                (*out) << qualityScore;

                // see if we mapped back to where we pulled the read from:
                if (se.getBestMatch().genomeMatchPosition == readGenomePosition)
                {
                    // simplest case - remapped back to the same location, so nothing else to print
                    (*out) << std::endl;   // all done with output
                }
                else
                {
                    std::string positionAndIndex;
                    // unless a repeat, we code it as 'M' -> mismapped
                    char reason = 'M';
                    // if valid map quality, but more than one postion, mark as a duplicate
                    if (se.getBestMatch().numBest > 1) reason = 'D';
                    gs->getChromosomeAndIndex(positionAndIndex, se.getBestMatch().genomeMatchPosition);
                    (*out) << "\t" << reason << "\t" << positionAndIndex << std::endl;;
                }
            }
            else
            {
                (*out) << "0\t";
                switch (se.getBestMatch().quality)
                {
                    case MatchedReadBase::UNSET_QUALITY:
                        (*out) << "Z";
                        break;
                    case MatchedReadBase::INVALID_DATA:
                        (*out) << "I";
                        break;
                    case MatchedReadBase::EARLYSTOP_QUALITY:
                        (*out) << "E";
                        break;
                    case MatchedReadBase::REPEAT_QUALITY:
                        (*out) << "R";
                        break;
                }
                (*out) << std::endl;
            }
        }
        (*out) << std::endl;
    }

    return 0;
}

int Test::test()
{
#if 0
    unsigned int word;
    uint64_t count = 0;

    // 10523 words were omitted due to exceeding 5000 cutoff.
    // 159356636 genome match positions were omitted due to those words exceeding the cutoff.
    for (word=0; word<wordIndex->wordsCountHashSize; word++)
    {
        if (wordIndex->wordReachedCutoff(word)) count++;
    }
    std::cout << count << " words reached cutoff." << std::endl;

    // this is pretty darn slow...
    wordIndex->testHashIndices(*gs);
#endif

    std::cout << "Karma unit tests: please add more, especially boundary cases!" << std::endl;

    testSEBaseSpaceReads();
    testPEBaseSpaceReads();

//    testRemapReference(30);     // 25, 35, 50, 75, 100
    return 0;
}

