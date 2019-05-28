#include "../tokenizer.hpp"
#include "gtest/gtest.h"

class TokenizerTest : public testing::Test {
protected:
    std::string line, option_name, option_contents;
    void catch_exception(std::string& lines);
};

void TokenizerTest::catch_exception(std::string& lines)
{
    try {
        TokenizeString(lines, option_name, option_contents);
        FAIL() << "Expected exception";
    }
    catch (const int& err) {
        EXPECT_EQ(err, -1);
    }
    catch (...) {
        FAIL() << "Expected exception to be int";
    }
}

TEST_F(TokenizerTest, ZeroContentStrings)
{
    line = "";
    EXPECT_FALSE(TokenizeString(line, option_name, option_contents));

    line = "  ";
    EXPECT_FALSE(TokenizeString(line, option_name, option_contents));

    line = " ()[]{}:,\t\n\n\v\f\r";
    EXPECT_FALSE(TokenizeString(line, option_name, option_contents));
}

TEST_F(TokenizerTest, CommentLine)
{
    line = "\%comment";
    EXPECT_FALSE(TokenizeString(line, option_name, option_contents));

    line = " \%comment";
    EXPECT_FALSE(TokenizeString(line, option_name, option_contents));

    line = " TEST = a \%end comment";
    EXPECT_TRUE(TokenizeString(line, option_name, option_contents));
}

TEST_F(TokenizerTest, ParsingIsCorrect)
{
    line = "TEST1=A1";
    EXPECT_TRUE(TokenizeString(line, option_name, option_contents));
    EXPECT_EQ(option_name, "TEST1");
    EXPECT_EQ(option_contents, "A1");

    line = "TEST2 = A2";
    EXPECT_TRUE(TokenizeString(line, option_name, option_contents));
    EXPECT_EQ(option_name, "TEST2");
    EXPECT_EQ(option_contents, "A2");

    line = "TEST3=A3 \%comment";
    EXPECT_TRUE(TokenizeString(line, option_name, option_contents));
    EXPECT_EQ(option_name, "TEST3");
    EXPECT_EQ(option_contents, "A3");

    line = " TEST4=A4";
    EXPECT_TRUE(TokenizeString(line, option_name, option_contents));
    EXPECT_EQ(option_name, "TEST4");
    EXPECT_EQ(option_contents, "A4");
}

TEST_F(TokenizerTest, Exceptions)
{
    line = "NO_EQUAL ";
    catch_exception(line);

    line = " = NO_ASSIGN_EQUAL ";
    catch_exception(line);

    line = "TWO NAMES = SOMETHING";
    catch_exception(line);

    line = "TWO = VALUES HERE";
    catch_exception(line);
}
