
#include "InplaceMerge.h"

#if defined(TEST)
#include "iostream"
#include "Generic.h"

int main(int argc, const char **argv)
{
    int a[] = {1,2,3,4,5};
    int b[] = {2,4,6,7,10};
    int c[] = {3,5,8,10,11};

    std::vector<int> z(15);

    std::copy(a, a+5, z.begin());
    std::copy(b, b+5, z.begin() + 5);
    std::copy(c, c+5, z.begin() + 10);

    std::vector<int> indeces, counts;

    indeces.push_back(0);
    indeces.push_back(5);
    indeces.push_back(10);

    counts.push_back(5);
    counts.push_back(5);
    counts.push_back(5);

    inplace_merge(indeces, counts, 0, 3, z);

    std::cout << z;
}

#endif
