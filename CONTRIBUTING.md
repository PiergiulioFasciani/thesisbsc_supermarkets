# Contributing Guide

## Development Setup

### Prerequisites
- R >= 4.0
- Python >= 3.8
- Git
- osmium CLI tool

### Installation Steps

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/thesis-codebase.git
   cd thesis-codebase
   ```

2. **Set up Python virtual environment:**
   ```bash
   python -m venv venv
   source venv/bin/activate  # macOS/Linux
   # or
   venv\Scripts\activate  # Windows
   ```

3. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

4. **Install R dependencies:**
   ```bash
   Rscript -e "install.packages(c('fixest', 'tidyverse', 'data.table', 'sf', 'yaml'))"
   Rscript scripts/install_honestdid.R
   ```

5. **Install system tools:**
   ```bash
   # macOS
   brew install osmium-tool
   
   # Ubuntu/Debian
   sudo apt-get install osmium-tool
   ```

## Code Style Guidelines

### R Code
- Use snake_case for variable names
- Use 2-space indentation
- Comment complex logic
- Use tidyverse/dplyr conventions for data manipulation

### Python Code
- Follow PEP 8 conventions
- Use 4-space indentation
- Add docstrings to functions
- Use type hints where practical

### Shell Scripts
- Use 2-space indentation
- Add error handling with `set -e`
- Document usage at script top

## Workflow for Modifications

1. **Create a feature branch:**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make changes** following code style guidelines

3. **Test locally:**
   ```bash
   # Run individual steps
   ./run_sampler.sh
   ./run_panel_builder.sh
   ./run_stacked_ppml_event_study.sh
   ```

4. **Update documentation** if behavior changes

5. **Commit with clear messages:**
   ```bash
   git commit -m "Add feature: clear description of change"
   ```

6. **Push and create Pull Request:**
   ```bash
   git push origin feature/your-feature-name
   ```

## Testing

### Validate R Scripts
```bash
Rscript -e "parse(file='scripts/stacked_ppml_event_study.R')"
```

### Validate Python Scripts
```bash
python -m py_compile scripts/osm_sampler.py scripts/panel_builder.py
```

### Test with Sample Data
```bash
# Use the provided Lombardy OSM file
./run_complete_analysis.sh
```

## Adding New Features

### Adding a New Analysis Type
1. Create new R script in `scripts/`
2. Add configuration section to `config.yml`
3. Create wrapper shell script `run_<feature>.sh`
4. Document in README.md
5. Update main orchestrator if needed

### Adding New Parameters
1. Update `config.yml` with new parameter and description
2. Update relevant analysis script to read from config
3. Update README.md Configuration section
4. Add validation for parameter values

### Adding New Outcomes
1. Update `outcome_types` in `config.yml`
2. Modify `panel_builder.py` OSM tag mapping if needed
3. Test panel generation with new outcomes
4. Results files will automatically include new outcomes

## Documentation

- Keep README.md synchronized with code changes
- Add comments to complex algorithms
- Document all configuration options
- Include example outputs in issues/PRs

## Reporting Issues

When reporting bugs, include:
- Steps to reproduce
- Expected vs. actual behavior
- `config.yml` settings used
- Output from running relevant scripts
- Error messages or logs

## Questions?

Open an issue or contact the maintainer.

---

Thank you for contributing!
