#!/usr/bin/env python3
import markdown

# Leggi il markdown
with open('POMPEJA_COMPARISON_REPORT_CORRECTED.md', 'r', encoding='utf-8') as f:
    md_content = f.read()

# Converti in HTML
html_body = markdown.markdown(md_content, extensions=['extra', 'codehilite', 'tables'])

# Template HTML completo
html = """<!DOCTYPE html>
<html lang="it">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pompeja Differential Correction Report</title>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            max-width: 900px;
            margin: 40px auto;
            padding: 20px;
            background: #f5f5f5;
        }
        h1 {
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }
        h2 {
            color: #34495e;
            margin-top: 30px;
            border-bottom: 2px solid #95a5a6;
            padding-bottom: 5px;
        }
        h3 {
            color: #555;
            margin-top: 20px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }
        th {
            background-color: #3498db;
            color: white;
        }
        tr:nth-child(even) {
            background-color: #f2f2f2;
        }
        code {
            background: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
        }
        pre {
            background: #2c3e50;
            color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
        }
        pre code {
            background: none;
            color: #ecf0f1;
        }
        .summary-box {
            background: #e8f5e9;
            border-left: 4px solid #4caf50;
            padding: 15px;
            margin: 20px 0;
        }
        .warning-box {
            background: #fff3e0;
            border-left: 4px solid #ff9800;
            padding: 15px;
            margin: 20px 0;
        }
        @media print {
            body {
                background: white;
                margin: 0;
                padding: 20px;
            }
            pre {
                page-break-inside: avoid;
            }
        }
    </style>
</head>
<body>
    """ + html_body + """
</body>
</html>"""

# Scrivi HTML
with open('POMPEJA_COMPARISON_REPORT_CORRECTED.html', 'w', encoding='utf-8') as f:
    f.write(html)

print('✓ Generated: POMPEJA_COMPARISON_REPORT_CORRECTED.html')
