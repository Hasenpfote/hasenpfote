#include "unicode.h"

namespace hasenpfote{

// utf8 to utf16
// gcc7.0.0 ... ---
// gcc6.1.0 ... ---
// gcc5.3.0 ... ---
// gcc5.2.0 ... need endianness.
// gcc5.1.0 ... ---
// clang4.0.0 . ---
std::u16string utf8_to_utf16(const std::string& s)
{
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
    // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
    std::wstring_convert<std::codecvt_utf8_utf16<std::uint16_t>, std::uint16_t> conv;
    auto temp = conv.from_bytes(s);
    return std::u16string(temp.cbegin(), temp.cend());
#else
    std::wstring_convert<std::codecvt_utf8_utf16<char16_t>, char16_t> conv;
    return conv.from_bytes(s);
#endif
}

// utf16 to utf8
// gcc7.0.0 ... ---
// gcc6.1.0 ... ---
// gcc5.3.0 ... ---
// gcc5.2.0 ... ---
// gcc5.1.0 ... ---
// clang4.0.0 . ---
std::string utf16_to_utf8(const std::u16string& s)
{
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
    // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
    std::wstring_convert<std::codecvt_utf8_utf16<std::uint16_t>, std::uint16_t> conv;
    auto p = reinterpret_cast<const std::uint16_t*>(s.c_str());
    return conv.to_bytes(p, p + s.length());
#else
    std::wstring_convert<std::codecvt_utf8_utf16<char16_t>, char16_t> conv;
    return conv.to_bytes(s);
#endif
}

// utf8 to utf32
// gcc7.0.0 ... ---
// gcc6.1.0 ... ---
// gcc5.3.0 ... ---
// gcc5.2.0 ... ---
// gcc5.1.0 ... ---
// clang4.0.0 . ---
std::u32string utf8_to_utf32(const std::string& s)
{
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
    // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
    std::wstring_convert<std::codecvt_utf8<std::uint32_t>, std::uint32_t> conv;
    auto temp = conv.from_bytes(s);
    return std::u32string(temp.cbegin(), temp.cend());
#else
    std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> conv;
    return conv.from_bytes(s);
#endif
}

// utf32 to utf8
// gcc7.0.0 ... ---
// gcc6.1.0 ... ---
// gcc5.3.0 ... ---
// gcc5.2.0 ... ---
// gcc5.1.0 ... ---
// clang4.0.0 . ---
std::string utf32_to_utf8(const std::u32string& s)
{
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
    // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
    std::wstring_convert<std::codecvt_utf8<std::uint32_t>, std::uint32_t> conv;
    auto p = reinterpret_cast<const std::uint32_t*>(s.c_str());
    return conv.to_bytes(p, p + s.length());
#else
    std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> conv;
    return conv.to_bytes(s);
#endif
}

// utf16 to utf32
// gcc7.0.0 ... need endianness.
// gcc6.1.0 ... need endianness.
// gcc5.3.0 ... need endianness.
// gcc5.2.0 ... need endianness.
// gcc5.1.0 ... need endianness.
// clang4.0.0 . need endianness.
std::u32string utf16_to_utf32(const std::u16string &s, bool is_big_endian)
{
    if(is_big_endian){
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
        // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
        std::wstring_convert<std::codecvt_utf16<std::uint32_t>, std::uint32_t> conv;
        const char16_t* data = s.c_str();
        auto bytes = conv.from_bytes(reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data + s.length()));
        return std::u32string(bytes.cbegin(), bytes.cend());
#else
        std::wstring_convert<std::codecvt_utf16<char32_t>, char32_t> conv;
        const char16_t* data = s.c_str();
        return conv.from_bytes(reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data + s.length()));
#endif
    }
    else{
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
        // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
        std::wstring_convert<std::codecvt_utf16<std::uint32_t, 0x10ffff, std::codecvt_mode::little_endian>, std::uint32_t> conv;
        const char16_t* data = s.c_str();
        auto bytes = conv.from_bytes(reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data + s.length()));
        return std::u32string(bytes.cbegin(), bytes.cend());
#else
        std::wstring_convert<std::codecvt_utf16<char32_t, 0x10ffff, std::codecvt_mode::little_endian>, char32_t> conv;
        const char16_t* data = s.c_str();
        return conv.from_bytes(reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data + s.length()));
#endif
    }
}

// utf32 to utf16
// gcc7.0.0 ... need endianness.
// gcc6.1.0 ... need endianness.
// gcc5.3.0 ... need endianness.
// gcc5.2.0 ... need endianness.
// gcc5.1.0 ... ---
// clang4.0.0 . need endianness.
std::u16string utf32_to_utf16(const std::u32string& s, bool is_big_endian)
{
    if(is_big_endian){
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
        // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
        std::wstring_convert<std::codecvt_utf16<std::uint32_t>, std::uint32_t> conv;
        auto p = reinterpret_cast<const std::uint32_t*>(s.c_str());
        auto bytes = conv.to_bytes(p, p + s.length());
        return std::u16string(reinterpret_cast<const char16_t*>(bytes.c_str()), bytes.length() / sizeof(char16_t));
#else
        std::wstring_convert<std::codecvt_utf16<char32_t>, char32_t> conv;
        auto bytes = conv.to_bytes(s);
        return std::u16string(reinterpret_cast<const char16_t*>(bytes.c_str()), bytes.length() / sizeof(char16_t));
#endif
    }
    else{
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
        // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
        std::wstring_convert<std::codecvt_utf16<std::uint32_t, 0x10ffff, std::codecvt_mode::little_endian>, std::uint32_t> conv;
        auto p = reinterpret_cast<const std::uint32_t*>(s.c_str());
        auto bytes = conv.to_bytes(p, p + s.length());
        return std::u16string(reinterpret_cast<const char16_t*>(bytes.c_str()), bytes.length() / sizeof(char16_t));
#else
        std::wstring_convert<std::codecvt_utf16<char32_t, 0x10ffff, std::codecvt_mode::little_endian>, char32_t> conv;
        auto bytes = conv.to_bytes(s);
        return std::u16string(reinterpret_cast<const char16_t*>(bytes.c_str()), bytes.length() / sizeof(char16_t));
#endif
    }
}

// ucs2 to utf8
// gcc7.0.0 ... ---
// gcc6.1.0 ... ---
// gcc5.3.0 ... ---
// gcc5.2.0 ... ---
// gcc5.1.0 ... ---
// clang4.0.0 . ---
std::string ucs2_to_utf8(const std::u16string& s)
{
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
    // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
    std::wstring_convert<std::codecvt_utf8<std::uint16_t>, std::uint16_t> conv;
    auto p = reinterpret_cast<const std::uint16_t*>(s.c_str());
    return conv.to_bytes(p, p + s.length());
#else
    std::wstring_convert<std::codecvt_utf8<char16_t>, char16_t> conv;
    return conv.to_bytes(s);
#endif
}

// utf8 to ucs2
// gcc7.0.0 ... ---
// gcc6.1.0 ... need endianness.
// gcc5.3.0 ... need endianness.
// gcc5.2.0 ... need endianness.
// gcc5.1.0 ... ---
// clang4.0.0 . ---
std::u16string utf8_to_ucs2(const std::string& s)
{
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
    // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
    std::wstring_convert<std::codecvt_utf8<std::uint16_t>, std::uint16_t> conv;
    auto temp = conv.from_bytes(s);
    return std::u16string(temp.cbegin(), temp.cend());
#else
    std::wstring_convert<std::codecvt_utf8<char16_t>, char16_t> conv;
    return conv.from_bytes(s);
#endif
}

// ucs2 to utf16
// gcc7.0.0 ... need endianness.
// gcc6.1.0 ... need endianness.
// gcc5.3.0 ... need endianness.
// gcc5.2.0 ... need endianness.
// gcc5.1.0 ... need endianness.
// clang4.0.0 . need endianness.
std::u16string ucs2_to_utf16(const std::u16string &s, bool is_big_endian)
{
    if(is_big_endian){
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
        // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
        std::wstring_convert<std::codecvt_utf16<std::uint16_t>, std::uint16_t> conv;
        auto p = reinterpret_cast<const std::uint16_t*>(s.c_str());
        auto bytes = conv.to_bytes(p, p + s.length());
        return std::u16string(reinterpret_cast<const char16_t*>(bytes.c_str()), bytes.length() / sizeof(char16_t));
#else
        std::wstring_convert<std::codecvt_utf16<char16_t>, char16_t> conv;
        auto bytes = conv.to_bytes(s);
        return std::u16string(reinterpret_cast<const char16_t*>(bytes.c_str()), bytes.length() / sizeof(char16_t));
#endif
    }
    else{
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
        // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
        std::wstring_convert<std::codecvt_utf16<std::uint16_t, 0x10ffff, std::codecvt_mode::little_endian>, std::uint16_t> conv;
        auto p = reinterpret_cast<const std::uint16_t*>(s.c_str());
        auto bytes = conv.to_bytes(p, p + s.length());
        return std::u16string(reinterpret_cast<const char16_t*>(bytes.c_str()), bytes.length() / sizeof(char16_t));
#else
        std::wstring_convert<std::codecvt_utf16<char16_t, 0x10ffff, std::codecvt_mode::little_endian>, char16_t> conv;
        auto bytes = conv.to_bytes(s);
        return std::u16string(reinterpret_cast<const char16_t*>(bytes.c_str()), bytes.length() / sizeof(char16_t));
#endif
    }
}

// utf16 to ucs2
// gcc7.0.0 ... need endianness.
// gcc6.1.0 ... need endianness.
// gcc5.3.0 ... need endianness.
// gcc5.2.0 ... need endianness.
// gcc5.1.0 ... need endianness.
// clang4.0.0 . need endianness.
std::u16string utf16_to_ucs2(const std::u16string &s, bool is_big_endian)
{
    if(is_big_endian){
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
        // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
        std::wstring_convert<std::codecvt_utf16<std::uint16_t>, std::uint16_t> conv;
        const char16_t* data = s.c_str();
        auto bytes = conv.from_bytes(reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data + s.length()));
        return std::u16string(bytes.cbegin(), bytes.cend());
#else
        std::wstring_convert<std::codecvt_utf16<char16_t>, char16_t> conv;
        const char16_t* data = s.c_str();
        auto temp = conv.from_bytes(reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data + s.length()));
        return std::u16string(temp.cbegin(), temp.cend());
#endif
    }
    else{
#if defined(_MSC_VER) && (_MSC_VER <= 1900)
        // Workaround for missing char16_t/char32_t instantiations in MSVC2015.
        std::wstring_convert<std::codecvt_utf16<std::uint16_t, 0x10ffff, std::codecvt_mode::little_endian>, std::uint16_t> conv;
        const char16_t* data = s.c_str();
        auto bytes = conv.from_bytes(reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data + s.length()));
        return std::u16string(bytes.cbegin(), bytes.cend());
#else
        std::wstring_convert<std::codecvt_utf16<char16_t, 0x10ffff, std::codecvt_mode::little_endian>, char16_t> conv;
        const char16_t* data = s.c_str();
        auto temp = conv.from_bytes(reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data + s.length()));
        return std::u16string(temp.cbegin(), temp.cend());
#endif
    }
}

}