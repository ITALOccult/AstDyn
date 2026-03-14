#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <iomanip>

int main() {
    std::ifstream f("/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp", std::ios::binary);
    if (!f) { std::cerr << "Fail open\n"; return 1; }
    
    char buf[1024];
    f.read(buf, 1024);
    
    int nd, ni, fwd, bwd;
    std::memcpy(&nd, buf + 8, 4);
    std::memcpy(&ni, buf + 12, 4);
    std::memcpy(&fwd, buf + 76, 4);
    std::memcpy(&bwd, buf + 80, 4);
    
    std::cout << "Header: ND=" << nd << " NI=" << ni << " FWD=" << fwd << " BWD=" << bwd << std::endl;
    
    if (fwd > 0) {
        f.seekg((fwd - 1) * 1024, std::ios::beg);
        f.read(buf, 1024);
        double next, prev, ns;
        std::memcpy(&next, buf, 8);
        std::memcpy(&prev, buf + 8, 8);
        std::memcpy(&ns, buf + 16, 8);
        std::cout << "Summary Rec " << fwd << ": Next=" << next << " Prev=" << prev << " NS=" << ns << std::endl;
        
        for (int i=0; i < (int)ns; ++i) {
            int offset = 24 + i * (nd * 8 + ni * 4);
            double t1, t2;
            std::memcpy(&t1, buf + offset, 8);
            std::memcpy(&t2, buf + offset + 8, 8);
            int body;
            std::memcpy(&body, buf + offset + 16, 4);
            std::cout << "  Seg " << i << ": Body=" << body << " T=[" << t1 << " : " << t2 << "]" << std::endl;
        }
        
        if (next > 0) {
            f.seekg((static_cast<int>(next) - 1) * 1024, std::ios::beg);
            f.read(buf, 1024);
            std::memcpy(&next, buf, 8);
            std::memcpy(&ns, buf + 16, 8);
            std::cout << "Summary Rec (Next): Next=" << next << " NS=" << ns << std::endl;
        }
    }
    
    return 0;
}
