#ifndef _MAIN_H_
#define _MAIN_H_

void mainCheck(const char *program, int argc, const char **argv);
void mainMap(const char *program, int argc, const char **argv);
void mainCreate(const char *program, int argc, const char **argv);
void mainHeader(const char *program, int argc, const char **argv);
void mainRemap(const char *program, int argc, const char **argv);
void mainTest(const char *program, int argc, const char **argv);
void commandList(std::ostream &out, int argc, const char **argv);
static std::string getReferenceNameWithArgs(
                                            int wordSize,
                                            int occurrenceCutoff,
                                            GenomeSequence *baseSpaceReference,
                                            GenomeSequence *colorSpaceReference);

#endif /* _MAIN_H_ */
