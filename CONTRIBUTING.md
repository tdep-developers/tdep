Contributing
===

When contributing to this repository, please first discuss the change you wish to make via issue, email, slack, or any other method with the maintainers of this repository. This will make life easier for everyone.

## Report Issues

Please use the [issue tracker](https://github.com/tdep-developers/tdep/issues) to report issues. Please try to answer these questions:

- Has this issue been discussed before? Please have a quick look at the existing issues. If not:
- What is the issue? What is the expected behavior?
- Is the problem reproducible? Please provide a _minimal_ example.


## Contribute Code via Pull Request

In order to contribute code to `TDEP`, please follow the usual steps for [preparing and creating a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) via your own [fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks). A few remarks regarding our guidelines for code and code style:

- Some of us use [`fprettify`](https://github.com/pseewald/fprettify) to format the code, the respective command is:
  ```bash
  fprettify -i 4 -l 500 --disable-indent-mod file.f90
  ```

  This will format the code with an indent (`-i`) of 4 characters, without enforcing line breaks (`-l 500` allows up to 500 characters per line - please use less.)

- Please _document_ and _test_ your changes by providing an example in the `examples` folder. 