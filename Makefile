.PHONY: gh-pages

GH_PAGES_BRANCH = gh-pages
MAIN_BRANCH = master
DOCS_DIR = docs
COMMIT_MSG = Update pkgdown site

gh-pages:
	@git rev-parse --is-inside-work-tree >/dev/null 2>&1 || (echo "Not a git repo" && exit 1)
	@test -d $(DOCS_DIR) || (echo "$(DOCS_DIR) not found" && exit 1)

	@echo "Switching to $(GH_PAGES_BRANCH)..."
	@git switch $(GH_PAGES_BRANCH) || git switch -c $(GH_PAGES_BRANCH)

	@echo "Cleaning branch..."
	@git rm -rf . >/dev/null 2>&1 || true

	@echo "Copying docs..."
	@cp -R $(DOCS_DIR)/* .

	@echo "Committing..."
	@git add .
	@git commit -m "$(COMMIT_MSG)" || echo "Nothing to commit"

	@echo "Pushing to origin/$(GH_PAGES_BRANCH)..."
	@git push -u origin $(GH_PAGES_BRANCH)

	@echo "Switching back to $(MAIN_BRANCH)..."
	@git switch $(MAIN_BRANCH)
