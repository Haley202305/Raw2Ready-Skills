# Raw2Ready Skills

A collection of **QoderWork Agent Skills** for bioinformatics workflows — from raw data to analysis-ready datasets.

Each skill is a self-contained module that teaches an AI agent how to perform a specific task, complete with utility scripts and step-by-step instructions.

## Available Skills

| Skill | Description |
|-------|-------------|
| [geo-matrix-merge](skills/geo-matrix-merge/) | Merge GEO supplementary matrix files (GSM\*.matrix.txt.gz) into a single AnnData (.h5ad) object for single-cell analysis. |

*More skills will be added over time.*

## How to Install a Skill

### Option 1 — Copy manually

Download the skill directory (e.g. `skills/geo-matrix-merge/`) and place it in your QoderWork skills folder:

| OS | Path |
|----|------|
| macOS / Linux | `~/.qoderworkcn/skills/geo-matrix-merge/` |
| Windows | `%USERPROFILE%\.qoderworkcn\skills\geo-matrix-merge\` |

### Option 2 — Git clone + symlink

```bash
git clone https://github.com/Haley202305/Raw2Ready-Skills.git
# Then symlink or copy individual skill directories to ~/.qoderworkcn/skills/
```

## Skill Structure

Every skill follows this layout:

```
skill-name/
├── SKILL.md          # Agent instructions (required)
└── scripts/          # Utility scripts (optional)
    └── tool.py
```

The agent reads `SKILL.md` to learn **what** to do and **when**, then calls scripts in `scripts/` to execute the work.

## Requirements

Individual skills list their own dependencies in their `SKILL.md`. Common dependencies across skills:

```bash
pip install pandas numpy anndata scipy
```

## Contributing

Pull requests welcome! To add a new skill:

1. Create a new directory under `skills/your-skill-name/`
2. Add a `SKILL.md` with YAML frontmatter (`name` + `description`)
3. Add utility scripts under `scripts/`
4. Update this README's skill table
5. Submit a PR

## License

MIT
