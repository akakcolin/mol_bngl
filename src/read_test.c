#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <regex.h>

// macro to determine nmatch from pmatch
#define ARRAY_SIZE(arr) (sizeof((arr)) / sizeof((arr)[0]))

// initialize target string and regular expression string
static const char *const str = "1) John Driverhacker;\n2) John Doe;\n3) John Foo;\n";
static const char *const re = "John.*o";

int main(void)
{
    // pointer to target string and regular expression
    static const char *s = str;
    regex_t     regex;

    // initialize pmatch containter for caching match offest and length
    regmatch_t  pmatch[1];

    // initialize offset and length types
    regoff_t    off, len;

    // exit if regular expression cant be compiled
    if (regcomp(&regex, re, REG_NEWLINE))
        exit(EXIT_FAILURE);

    printf("String = \"%s\"\n", str);
    printf("Matches:\n");

    // for each possible match
    for (int i = 0; ; i++) {

        // exit if no more matches
        if (regexec(&regex, s, ARRAY_SIZE(pmatch), pmatch, 0))
            break;

        // compute offset of match and length of match and print
        off = pmatch[0].rm_so + (s - str);
        len = pmatch[0].rm_eo - pmatch[0].rm_so;
        printf("#%d:\n", i);
        printf("offset = %jd; length = %jd\n", (intmax_t) off, (intmax_t) len);

        // print the match
        printf("substring = \"%.*s\"\n", len, s + pmatch[0].rm_so);

        // move the pointer to the next start of the string
        s += pmatch[0].rm_eo;
    }

    exit(EXIT_SUCCESS);
}#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <regex.h>

// macro to determine nmatch from pmatch
#define ARRAY_SIZE(arr) (sizeof((arr)) / sizeof((arr)[0]))

// initialize target string and regular expression string
static const char *const str = "ISOPS";
static const char *const re = "([A-Z]+)";

int main(void)
{
    // pointer to target string and regular expression
    static const char *s = str;
    regex_t     regex;

    // initialize pmatch containter for caching match offest and length
    regmatch_t  pmatch[1];

    // initialize offset and length types
    regoff_t    off, len;

    // exit if regular expression cant be compiled
    if (regcomp(&regex, re, REG_NEWLINE))
        exit(EXIT_FAILURE);

    printf("String = \"%s\"\n", str);

    // for each possible match
    for (int i = 0; ; i++) {
        // exit if no more matches
        if (regexec(&regex, s, ARRAY_SIZE(pmatch), pmatch, 0))
            break;

        // compute offset of match and length of match and print
        off = pmatch[0].rm_so + (s - str);
        len = pmatch[0].rm_eo - pmatch[0].rm_so;
        printf("#%d:\n", i);
        printf("offset = %jd; length = %jd\n", (intmax_t) off, (intmax_t) len);

        // print the match
        printf("substring = \"%.*s\"\n", len, s + pmatch[0].rm_so);

        // move the pointer to the next start of the string
        s += pmatch[0].rm_eo;
    }

    exit(EXIT_SUCCESS);
}
