# Contribution Guidelines

So you are wondering how can you contribute to TARDIS? Congrats, you've landed on the right page!

## Imposter Syndrome Disclaimer

_We want your help. No really, we do._

There might be a little voice inside that tells you you're not ready; that you need to do one more tutorial, or learn another framework, or write a few more blog posts before you can help me with this project.

We assure you, that's not the case.

TARDIS has some clear "Contribution Guidelines" that you can read below in the following sections.

The contribution guidelines outline the process that you'll need to follow to get a patch merged. By making expectations and process explicit, we hope it will make it easier for you to contribute.

And you don't just have to write code. You can help out by writing documentation, tests, or even by giving feedback about this work. (And yes, that includes giving feedback about the contribution guidelines.)

(Thanks to [Adrienne Lowe](https://github.com/adriennefriend/imposter-syndrome-disclaimer) who came up with this disclaimer and yeah, made it OPEN for anyone to use!)

## How can I contribute?

There are multiple ways in which you can help us:

- Found a bug in TARDIS? Report it to us!
- Caught a typo in documentation or want to make it better to understand? Edit it!
- Know how to fix an issue or add a new feature? Make a patch!
- Loved using TARDIS? Share it with others!
- Anything we missed to mention? Then, what are you waiting for!

### Reporting a Bug

TARDIS is in active development. There's no surprise that you may encounter something that doesn't work for your use case. Or maybe you have some suggestions about how can we improve some functionality. Feel free to share any of it with us by [opening an issue](https://docs.github.com/en/github/managing-your-work-on-github/creating-an-issue) [here](https://github.com/tardis-sn/tardis/issues/).

Please make sure that you provide all the necessary information requested by prompts in the issue body - it will not only make our work easier but will also help you to communicate your problem better.

### Editing the Documentation

There is always a scope of improvement in documentation to add some missing information or to make it easier for reading. And here lies an opportunity for you, since you can edit the documentation page you want which is stored as a text file in [`docs`](https://github.com/tardis-sn/tardis/tree/master/docs) directory of TARDIS.

After editing the file locally, build the docs as described in [these instructions](https://tardis-sn.github.io/tardis/development/documentation_guidelines.html#building-documentation-locally) and then you can submit your changes to us by making a patch as described in the next section.

### Making a Patch

If you have peeked in our codebase and realized how you can fix a problem or if you know how to add a new feature, well done! If not, don't worry - just pick an [easy](https://github.com/tardis-sn/tardis/labels/easy) or [good-first](https://github.com/tardis-sn/tardis/labels/good%20first%20issue) issue and get started to fix it.

To contribute your code to TARDIS, you'll need to make a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) from your fork of TARDIS repository. This development workflow using Git may look daunting at first, but it is not if you follow [this guide](https://tardis-sn.github.io/tardis/development/git_workflow.html#preparation-and-working-with-git) that we have prepared for you.

When you make a pull request, please provide all the necessary information requested by prompts in the pull request body. Also, make sure that the code you're submitting always accounts for the following three:

- **Maintaining code quality:** Your code must follow the PEP8 style guide, should cover edge cases, etc. Check our [code quality guidelines](https://tardis-sn.github.io/tardis/development/code_quality.html) for more details.
- **Documenting the code:** You must write docstrings in functions/classes, put a relevant example in TARDIS docs and make sure docs get built correctly. This is explained in detail in our [documentation guidelines](https://tardis-sn.github.io/tardis/development/documentation_guidelines.html).
- **Testing the code:** There should be unit-tests for most of the functions/methods and they must pass our testing framework. To run test locally, you can follow [this guide](https://tardis-sn.github.io/tardis/development/running_tests.html).

### Spreading the word of mouth

If you find TARDIS helpful, you can share it with your peers, colleagues, and anyone who can benefit from TARDIS. If you've used TARDIS in your research, please make sure to cite us (https://tardis-sn.github.io/tardis/credits.html). By telling other people about how we helped you, you'll help us in turn, in extending our impact. And we would absolutely love it if you give us a shout-out on [twitter](https://twitter.com/tardis_sn/)!

## What if I need help?

Contact us at [Gitter](https://gitter.im/tardis-sn/tardis) without any hesitation. We already appreciate that you're trying to help us, so feel free to ask any doubt you have or whatever is keeping you stuck.

---

Thank you for contributing (because we believe in YOU)!
