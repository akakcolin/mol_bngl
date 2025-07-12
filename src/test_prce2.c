#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PCRE2_CODE_UNIT_WIDTH 8
#include <pcre2.h>

int main() {
    char *subject = "cooh!4,cooh!5,cooh!6";
    char *pattern = "([^, ]+)";
    int count = 0;
    int rc;
    int ovector[30];
    pcre2_code *re;
    pcre2_match_data *match_data;
    PCRE2_SIZE subject_length;
    PCRE2_SPTR pa=(PCRE2_SPTR)pattern;
    int match_count;
    PCRE2_SIZE erroffset;
    int errnum;

    re = pcre2_compile((PCRE2_SPTR)pattern, PCRE2_ZERO_TERMINATED, 0, &errnum, &erroffset, NULL);
    if (!re) {
        printf("PCRE2 compilation failed\n");
        return 1;
    }

    match_data = pcre2_match_data_create_from_pattern(re, NULL);
    if (!match_data) {
        printf("PCRE2 match data creation failed\n");
        pcre2_code_free(re);
        return 1;
    }

    subject_length = strlen(subject);
    printf("%d\n", subject_length);

    rc = pcre2_match(re, (PCRE2_SPTR)subject, subject_length, 0, 0, match_data, NULL);
    if (rc < 0) {
        printf("PCRE2 match failed\n");
        pcre2_match_data_free(match_data);
        pcre2_code_free(re);
        return 1;
    }

    match_count = pcre2_get_ovector_pointer(match_data)[0];
    printf("Found %d %d matches:\n", match_count, rc);


    while (rc > 0) {
        PCRE2_SIZE start_offset = pcre2_get_ovector_pointer(match_data)[0];
        PCRE2_SIZE end_offset = pcre2_get_ovector_pointer(match_data)[1];
        const char *substring = subject + start_offset;
        size_t substring_length = end_offset - start_offset;

        printf("Match %d: %.*s\n", count + 1, (int)substring_length, substring);

        rc = pcre2_match(re, (PCRE2_SPTR8)subject, subject_length, end_offset, 0, match_data, NULL);
        count++;
    }

    pcre2_match_data_free(match_data);
    pcre2_code_free(re);

    printf("Total number of matches: %d\n", count);

    return 0;
}
