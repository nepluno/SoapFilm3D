#pragma once

/* @class      : HAS_TYPEDEF
 * @brief      : This macro may be used to check if a class has a public typedef
 * @param NAME        : Name of struct this macro creates.
 * @param TYPEDEF     : Name of typedef to test for.
 * @note This macro is C++03 compatible
 */
#define HAS_TYPEDEF(NAME, TYPEDEF)                                      \
  template <typename CLASS>                                             \
  struct NAME {                                                         \
    template <typename U> static char chk(U::TYPEDEF*);                 \
    template <typename  > static long chk(...);                         \
    static const bool value = sizeof(chk<CLASS>(0)) == sizeof(char);    \
  }

/* @class      : HAS_MEM_FUNC
 * @brief      : This macro may be used to check if a class has a public
 *               const member function with particular signature.
 * @param NAME        : Name of struct this macro creates.
 * @param RETURN_TYPE : Return type of the member function to test for.
 * @param FUNC        : Name of member function to test for.
 * @param (...)       : The argument types of the member function to test for.
 *                      These complete the signature of the member funtion.
 * @note This macro is C++03 compatible
 */
#define HAS_MEM_FUNC(NAME, RETURN_TYPE, FUNC, ...)                      \
  template <typename CLASS>                                             \
  struct NAME {                                                         \
    typedef RETURN_TYPE (CLASS::*A)(__VA_ARGS__) const;                 \
    template <typename U, U> struct type_check;                         \
    template <typename U> static char chk(type_check<A, &U::FUNC>*);    \
    template <typename  > static long chk(...);                         \
    static bool const value = sizeof(chk<CLASS>(0)) == sizeof(char);    \
  }
