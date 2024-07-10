---
name: Error
about: Report a TDEP error message and ask for help
title: "[ERROR] Your error message"
labels: error
assignees: ''

---

Before reporting an error, please check if your problem was already covered in on of the following resources:

- [ ] [Troubleshooting section](https://github.com/tdep-developers/tdep?tab=readme-ov-file#troubleshooting)
- [ ] [Issue tracker](https://github.com/tdep-developers/tdep/issues) including [closed issues](https://github.com/tdep-developers/tdep/issues?q=is%3Aissue+is%3Aclosed)

If this did not help, please report your error detailed as possible.

Therefore, please provide the following:

- [ ] Report the `git commit` or release version of TDEP that you are working with.
- [ ] Describe the error message you encountered.
- [ ] Provide the **full command** you were running (e.g. `extract_forceconstants -rc2 5`)
  - [ ] If you are using a script to run or submit the command, please provide this script as well
- [ ] Add a **full** log file for the binary you were running (e.g. via ``extract_forceconstants -rc2 5 2>&1 | tee extract_forceconstants.log`)
- [ ] Add **all input files** necessary to run the command
