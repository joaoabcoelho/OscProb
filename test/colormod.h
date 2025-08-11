//#include <ostream>

using namespace std;

namespace Color {

    enum Code {
        FG_RED      = 31,
        FG_GREEN    = 32,
        FG_BLUE     = 34,
        FG_DEFAULT  = 39,
        BG_RED      = 41,
        BG_GREEN    = 42,
        BG_BLUE     = 44,
        BG_DEFAULT  = 49
    };

    string Modifier(string s, Code code) {
      return "\033[" + to_string(code) + "m" + s + "\033[0m";
    }

    static const string FAILED = Modifier("FAILED:", FG_RED);
    static const string PASSED = Modifier("PASSED:", FG_GREEN);

}
