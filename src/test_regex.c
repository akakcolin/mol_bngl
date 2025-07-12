// Example of using regex to parse "ISO(cooh,cooh)"
#include <regex.h>
#include <stdio.h>
/*
re_PG_part = re.compile(r'PG\(oh!*(\%?\d*),oh!*(\%?\d*)\)')
re_TMP_part = re.compile(r'TMP\(oh!*(\%?\d*),oh!*(\%?\d*),oh!*(\%?\d*)\)')
re_ISOPS_part = re.compile(r'ISOPS\(cooh!*(\%?\d*),cooh!*(\%?\d*)\)')



_ISOPS_smile='O=C(O[*])c1cccc(C(=O)O[*])c1'
_PG_smile = 'CC([*])C[*]'
_TMP_smile = 'CCC(C[*])(C[*])C[*]'
_H2O_smile = 'O'

_ISOPS_func_index0=[7,24]
_TMP_func_index0=[7,13,18]
_PG_func_index0=[5,10]
*/

char* display_match_re(const char* s, regmatch_t matches[], int mno) {
  static char buf[1024];
  snprintf(buf, sizeof buf, "%.*s", (int)(matches[mno].rm_eo - matches[mno].rm_so), s+matches[mno].rm_so);
  return buf;
}

int main(int argc, const char **argv) {
  regex_t regex, regex_comp;
  regmatch_t matches[3];  // Array to store match positions
  const char *pattern =
      "([A-Za-z]+)\\(([^)]+)\\)";  // Regex pattern to match "ISO(cooh,cooh)"
  const char *componentPattern="([^ ,]+)";
  const char *input = "ISO(cooh,cooh)";
  regmatch_t matches_comp[3];

  char comp[32];
  // Compile the regex pattern
  if (regcomp(&regex, pattern, REG_EXTENDED) != 0) {
    fprintf(stderr, "Failed to compile regex\n");
    return 0;
  }

  // Compile the regex pattern
  if (regcomp(&regex_comp, componentPattern , REG_EXTENDED) != 0) {
    fprintf(stderr, "Failed to compile regex comp\n");
    return 0;
  }

  // Execute the regex
  if (regexec(&regex, input, 3, matches, 0) == 0) {
    char name[32], groups[32];
    // Extract the matched groups
    snprintf(name, matches[1].rm_eo - matches[1].rm_so + 1, "%s",
             input + matches[1].rm_so);
    snprintf(groups, matches[2].rm_eo - matches[2].rm_so + 1, "%s",
             input + matches[2].rm_so);
    printf("Parsed name: %s, groups: %s\n", name, groups);

    if(regexec(&regex_comp, groups, 3, matches_comp, 0) == 0){
      snprintf(comp, matches_comp[1].rm_eo - matches_comp[1].rm_so + 1, "%s", groups + matches_comp[1].rm_so);
      printf("group %s\n", comp);

      snprintf(comp, matches_comp[1].rm_eo - matches_comp[1].rm_so + 1, "%s", groups + matches_comp[1].rm_so);

    }


  } else {
    fprintf(stderr, "No match found\n");
  }

  // Free the compiled regex
  regfree(&regex);
  return 0;
}
