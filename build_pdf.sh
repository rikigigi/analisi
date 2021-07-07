pandoc README_.md -V geometry:a4paper -o README.pdf
python -m readme2tex --output README.md README_.md  --nocdn
