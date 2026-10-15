# T4 Checkpoint 1 - Data start pack

- What this is: the 'start-evidence package' for the data part of Checkpoint 1 (download log + source, gene/label and split dictionaries).
- Yellow cells marked ENTER = you fill these in while downloading (release version, date, size, status, screenshot name).
- Step 1: open each URL in 'Download_Log', note the release/version shown on the page, start the download, screenshot it.
- Step 2: fill in date, size and status (Done / In progress / Not started). Unfinished big downloads are OK at Checkpoint 1.
- Step 3: the three dictionary sheets are drafts of what the tables WILL contain. Check them as a group; no data has to be computed yet.
- Key rule (leakage): nothing derived from DepMap gene effects may be used as a feature. DepMap only makes the label.
- Key rule (scope): the TROP2 paper reproduction (E-MTAB-16585) and the LUAD dependency model are two separate tasks.
- Release versions were not verified by Claude (sites were not reachable from its sandbox). Record what each site shows on the day you download.

## Where things live
- GitHub (`data/dictionaries/`): these CSVs and the code. Small text files only.
- Box: the downloaded data files (too large for GitHub) and download screenshots.
- Add `data/raw/` to `.gitignore` so big files never get committed.
