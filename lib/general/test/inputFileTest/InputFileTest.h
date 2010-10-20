#include <string>
#include "InputFile.h"


// This test class is only for when zlib is available.
#ifdef __ZLIB_AVAILABLE__

class IFILE_Test : public InputFile
{
public:
    void test();

    static const int TEST_FILE_SIZE;
    static const int BGZF_TEST_FILE_SIZE;
    static const std::string TEST_FILE_CONTENTS;

private:
    void test_readFromFile(const char* extension);

    // Tested together because they are used to test each other.
    void test_ifeof_ifrewind(const char* extension);

    // Tested together to verify they can be successfully be called after the
    // other has been called.
    void test_ifread_ifgetc(const char* extension);

    void test_ifclose(const char* extension);

    void test_ifseek(const char* extension);

    void test_noExistRead(const char *extension);

    void openFile(const char* extension);
    void openLargeFile(const char* extension);
    void openNoExistFile(const char* extension);

    // Buffer used for reading into.
    static const int MAX_TEST_BUFFER_SIZE = 100;
    char myTestBuffer[MAX_TEST_BUFFER_SIZE];

};

#endif
