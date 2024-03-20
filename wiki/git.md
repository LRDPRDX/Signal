# Git development workflow

## The `main` branch

This branch represents the most stable working version of the repository.

:warning: You should never commit directly to `main`. :warning:

## Creating a new feature branch

When you want to start working on a bug or feature you should (read _MUST_)
create a new branch --- _a feature branch_ starting from the latest revision
of `main` :

 - Create an issue in the GitHub interface with a reasonable title
 - Check it locally

When you're done with (a part of) your work commit your changes :

```
git add .
git commit
```

and push it to the remote :

```
git push
```

