import os


project = 'HarP'
copyright = '2026, Mattys Pouyez'
author = 'Mattys Pouyez'
version = "1.0"
release = "1.0.0"

extensions = []



templates_path = ['_templates']
exclude_patterns = []


numfig = True

html_theme_path = ['.']

html_theme = "smilei_theme"
html_logo  = "../logo/HarPLogoWhite.svg"

html_static_path = ['_static']

html_context = {
    "smilei_theme_version": "5.1",
    "logo": "HarPLogoWhite.svg"
}


