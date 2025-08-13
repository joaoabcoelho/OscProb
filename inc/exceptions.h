#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

//.............................................................................
///
/// @brief Helper function to split a string by a delimiter.
///
/// This is used to parse the comma-separated argument names.
///
/// @param s         - The string to split.
/// @param delimiter - The character delimiter.
///
/// @return A vector of the split strings.
///
inline std::vector<std::string> split(const std::string& s, char delimiter)
{
  std::string              token;
  std::vector<std::string> tokens;
  std::istringstream       tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    // Trim leading/trailing whitespace
    size_t first = token.find_first_not_of(" \t\n\r");
    if (std::string::npos != first) {
      size_t last = token.find_last_not_of(" \t\n\r");
      tokens.push_back(token.substr(first, (last - first + 1)));
    }
    else {
      tokens.push_back("");
    }
  }
  return tokens;
}

//.............................................................................
///
/// @brief The main variadic template function to log the arguments.
///
/// It unpacks the names and values from the macro call
/// and returns a string of name = value pairs.
///
/// @tparam Args The types of the arguments.
///
/// @param names - The comma-separated string of argument names.
/// @param args  - The actual argument values.
///
/// @return The string containing the formatted log message.
///
template <typename... Args>
std::string format_args(const std::string& names, const Args&... args)
{
  std::stringstream ss;

  std::vector<std::string> arg_names = split(names, ',');

  std::size_t i = 0;
  // Lambda to iterate and print each argument.
  auto printer = [&](const auto& arg_value) {
    if (i < arg_names.size()) {
      if (i > 0) ss << ", ";
      ss << arg_names[i] << " = " << arg_value;
    }
    i++;
  };
  // Unpack the arguments and print them using the lambda.
  (printer(args), ...);

  return ss.str();
}

// Macro to get the most descriptive function name available
#if defined(__GNUC__) || defined(__clang__)
#define FUNCTION_NAME __PRETTY_FUNCTION__
#elif defined(_MSC_VER)
#define FUNCTION_NAME __FUNCSIG__
#else
#define FUNCTION_NAME __func__
#endif

// The core macro for throwing an exception with detailed info
#define THROW_ON_INVALID_ARG(condition, ...)                                   \
  do {                                                                         \
    if (!(condition)) {                                                        \
      std::stringstream ss;                                                    \
      ss << "\n    Condition '" << #condition << "' failed."                   \
         << "\n    With: " << format_args(#__VA_ARGS__, __VA_ARGS__)           \
         << "\n    In function: " << FUNCTION_NAME                             \
         << "\n    At file: " << __FILE__ << "\n    At line: " << __LINE__;    \
      throw std::invalid_argument(ss.str());                                   \
    }                                                                          \
  }                                                                            \
  while (0)

#endif // EXCEPTIONS_H
