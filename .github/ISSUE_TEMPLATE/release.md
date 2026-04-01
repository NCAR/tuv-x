---
name: Version Release
about: Create an issue to make a new release
title: 'Release X.X.X'
labels: ''
assignees: ''

---

## Testing

- [ ] GitHub Actions are passing on `main`

## Deployment

- [ ] Create a new branch (do **not** name it `release`)
- [ ] Update the version number in `CMakeLists.txt`
- [ ] Update `CITATION.cff`
  - Update the version number
  - Ensure all contributors are listed as authors
- [ ] On GitHub, merge `main` into `release` — **do NOT squash and merge**
  - Alternatively, merge locally and push: `git checkout release && git merge main && git push`
- [ ] Create a tag and add release notes on GitHub
  - Be sure to select the `release` branch as the target
