#include <iostream>
#include "ShakeAVX2.h"

using namespace std;

int main()
{
    uint64_t nonce = 0x123456789abcdef;
    uint64_t counter = 1;

    cout << "Init" << endl;
    ShakeAVX2 shake(nonce, counter);
    cout << "Print states" << endl;
    shake.print_state();
    cout << "Update" << endl;
    shake.update(nonce, counter + 1);
    cout << "Print states" << endl;
    shake.print_state();

    return 0;
}
