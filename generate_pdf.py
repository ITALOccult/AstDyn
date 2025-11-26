#!/usr/bin/env python3
"""
Convert Markdown report to PDF using built-in tools
"""

import subprocess
import sys
import os

def md_to_pdf_mac(md_file, pdf_file):
    """Convert MD to PDF on macOS using textutil and cupsfilter"""
    
    # Step 1: MD → HTML (using Python markdown or simple conversion)
    html_file = md_file.replace('.md', '.html')
    
    with open(md_file, 'r', encoding='utf-8') as f:
        md_content = f.read()
    
    # Simple HTML wrapper with better styling
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            line-height: 1.6;
            max-width: 210mm;
            margin: 0 auto;
            padding: 20mm;
            font-size: 11pt;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
            page-break-before: always;
        }}
        h1:first-of-type {{
            page-break-before: avoid;
        }}
        h2 {{
            color: #34495e;
            border-bottom: 2px solid #95a5a6;
            padding-bottom: 5px;
            margin-top: 30px;
        }}
        h3 {{
            color: #7f8c8d;
            margin-top: 20px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 15px 0;
            font-size: 10pt;
        }}
        th, td {{
            border: 1px solid #bdc3c7;
            padding: 8px;
            text-align: left;
        }}
        th {{
            background-color: #ecf0f1;
            font-weight: bold;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        code {{
            background-color: #f4f4f4;
            padding: 2px 5px;
            border-radius: 3px;
            font-family: 'Monaco', 'Courier New', monospace;
            font-size: 9pt;
        }}
        pre {{
            background-color: #2c3e50;
            color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
            font-size: 9pt;
        }}
        pre code {{
            background-color: transparent;
            color: inherit;
        }}
        blockquote {{
            border-left: 4px solid #3498db;
            padding-left: 15px;
            color: #7f8c8d;
            margin: 15px 0;
        }}
        .check {{
            color: #27ae60;
            font-weight: bold;
        }}
        .cross {{
            color: #e74c3c;
            font-weight: bold;
        }}
        hr {{
            border: none;
            border-top: 2px solid #bdc3c7;
            margin: 30px 0;
        }}
        @media print {{
            body {{
                font-size: 10pt;
            }}
            h1 {{
                page-break-after: avoid;
            }}
            table {{
                page-break-inside: avoid;
            }}
        }}
    </style>
</head>
<body>
"""
    
    # Simple MD to HTML conversion (basic)
    lines = md_content.split('\n')
    in_code_block = False
    in_table = False
    
    for line in lines:
        # Code blocks
        if line.strip().startswith('```'):
            if in_code_block:
                html_content += '</code></pre>\n'
                in_code_block = False
            else:
                html_content += '<pre><code>'
                in_code_block = True
            continue
        
        if in_code_block:
            html_content += line + '\n'
            continue
        
        # Headers
        if line.startswith('# '):
            html_content += f'<h1>{line[2:]}</h1>\n'
        elif line.startswith('## '):
            html_content += f'<h2>{line[3:]}</h2>\n'
        elif line.startswith('### '):
            html_content += f'<h3>{line[4:]}</h3>\n'
        elif line.startswith('#### '):
            html_content += f'<h4>{line[5:]}</h4>\n'
        
        # Horizontal rule
        elif line.strip() == '---':
            html_content += '<hr>\n'
        
        # Tables (simple detection)
        elif '|' in line and not line.startswith('<!'):
            if not in_table:
                html_content += '<table>\n'
                in_table = True
            
            cells = [c.strip() for c in line.split('|')[1:-1]]
            
            # Header separator line
            if all(c.replace('-', '').strip() == '' for c in cells if c):
                continue
            
            # Detect if it's header row (first row of table)
            is_header = '---' in '\n'.join(lines[lines.index(line):lines.index(line)+2])
            
            html_content += '<tr>\n'
            tag = 'th' if is_header else 'td'
            for cell in cells:
                # Replace checkmarks
                cell = cell.replace('✅', '<span class="check">✓</span>')
                cell = cell.replace('❌', '<span class="cross">✗</span>')
                cell = cell.replace('⭐', '★')
                html_content += f'<{tag}>{cell}</{tag}>\n'
            html_content += '</tr>\n'
        
        else:
            if in_table and '|' not in line:
                html_content += '</table>\n'
                in_table = False
            
            # Regular paragraphs
            if line.strip():
                # Replace markdown bold/italic
                line = line.replace('**', '<strong>').replace('**', '</strong>')
                line = line.replace('*', '<em>').replace('*', '</em>')
                # Replace checkmarks
                line = line.replace('✅', '<span class="check">✓</span>')
                line = line.replace('❌', '<span class="cross">✗</span>')
                line = line.replace('⚠️', '⚠')
                line = line.replace('⭐', '★')
                
                # Inline code
                import re
                line = re.sub(r'`([^`]+)`', r'<code>\1</code>', line)
                
                html_content += f'<p>{line}</p>\n'
    
    if in_table:
        html_content += '</table>\n'
    
    html_content += """
</body>
</html>
"""
    
    # Write HTML
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✓ Generated HTML: {html_file}")
    
    # Step 2: HTML → PDF using macOS built-in tools
    try:
        # Try cupsfilter (best quality)
        subprocess.run([
            'cupsfilter',
            html_file,
            '>', pdf_file
        ], check=True, shell=True)
        print(f"✓ Generated PDF: {pdf_file}")
        return True
    except:
        pass
    
    try:
        # Try wkhtmltopdf if available
        subprocess.run([
            'wkhtmltopdf',
            '--enable-local-file-access',
            '--page-size', 'A4',
            '--margin-top', '20mm',
            '--margin-bottom', '20mm',
            html_file,
            pdf_file
        ], check=True)
        print(f"✓ Generated PDF: {pdf_file}")
        return True
    except:
        pass
    
    # Fallback: Just keep HTML
    print(f"⚠ PDF generation failed, but HTML is available: {html_file}")
    print(f"  You can open the HTML file in Safari and print to PDF manually")
    return False

if __name__ == '__main__':
    md_file = 'POMPEJA_COMPARISON_REPORT.md'
    pdf_file = 'POMPEJA_COMPARISON_REPORT.pdf'
    
    if not os.path.exists(md_file):
        print(f"Error: {md_file} not found")
        sys.exit(1)
    
    success = md_to_pdf_mac(md_file, pdf_file)
    
    if success:
        print("\n✅ PDF report generated successfully!")
    else:
        print("\n⚠️  HTML report generated, manual PDF conversion needed")
    
    sys.exit(0 if success else 1)
