# Contributing

Thank you for your interest in contributing to MACE Inference!

## Ways to Contribute

- 🐛 Report bugs
- 💡 Suggest features
- 📝 Improve documentation
- 🔧 Submit pull requests

## Development Setup

1. **Fork and clone the repository**

    ```bash
    git clone https://github.com/YOUR_USERNAME/mace-inference.git
    cd mace-inference
    ```

2. **Create a virtual environment**

    ```bash
    python -m venv venv
    source venv/bin/activate  # Linux/macOS
    # or: venv\Scripts\activate  # Windows
    ```

3. **Install in development mode**

    ```bash
    pip install -e ".[dev]"
    ```

4. **Run tests to verify setup**

    ```bash
    pytest tests/ -v
    ```

## Code Style

We use the following tools:

- **Ruff** for linting
- **Black** for formatting (optional)
- **mypy** for type checking

Run linting:

```bash
ruff check src/
```

Run type checking:

```bash
mypy src/mace_inference --ignore-missing-imports
```

## Testing

Run the test suite:

```bash
# All tests
pytest tests/ -v

# Specific test file
pytest tests/test_core.py -v

# With coverage
pytest tests/ --cov=mace_inference --cov-report=html
```

## Pull Request Process

1. **Create a feature branch**

    ```bash
    git checkout -b feature/your-feature-name
    ```

2. **Make your changes**
    - Write clear, documented code
    - Add tests for new features
    - Update documentation if needed

3. **Run tests and linting**

    ```bash
    pytest tests/ -v
    ruff check src/
    ```

4. **Commit with clear messages**

    ```bash
    git commit -m "Add feature: description of what you added"
    ```

5. **Push and create PR**

    ```bash
    git push origin feature/your-feature-name
    ```

    Then create a Pull Request on GitHub.

## Commit Message Guidelines

Use clear, descriptive commit messages:

- `Add feature: bulk modulus calculation`
- `Fix: handle empty atoms object in single_point`
- `Docs: update installation guide for Windows`
- `Test: add tests for D3 correction`
- `Refactor: simplify device detection logic`

## Adding New Features

When adding a new calculation type:

1. **Add the task function** in `src/mace_inference/tasks/`
2. **Add the method** to `MACEInference` in `core.py`
3. **Add tests** in `tests/`
4. **Add documentation** in `docs/`
5. **Add an example** in `examples/`

## Documentation

Documentation is built with MkDocs. To preview locally:

```bash
pip install mkdocs-material mkdocstrings[python]
mkdocs serve
```

Then open http://localhost:8000

## Reporting Issues

When reporting bugs, please include:

1. Python version and OS
2. MACE Inference version
3. Minimal code to reproduce
4. Full error traceback
5. What you expected vs. what happened

## Code of Conduct

- Be respectful and inclusive
- Focus on constructive feedback
- Help others learn and grow

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

## Questions?

Open an issue on GitHub or start a discussion!
