#!/usr/bin/env python3
import os
import re
import sys

# Rules from "Code Fits in Your Head"
RULES = {
    "R1": "Function exceeds 15-20 lines (Method should fit on one screen)",
    "R2": "Variable name suggests implementation instead of physical meaning",
    "R3": "Comment repeating code logic",
    "R4": "Potential breach of Single Responsibility Principle (Multi-tasking)",
    "R5": "Magic number detected",
    "R6": "Raw 'double' used instead of Physical Type (Distance, Velocity, etc.)",
    "R7": "Silent failure or empty pointer check detected",
    "R8": "Too many boolean flags in function parameters",
}

# Regex Patterns
RE_FUNCTION = re.compile(r'(?:[\w:<>&]+\s+)+(\w+)\s*\(([^)]*)\)\s*(?:const\s*)?\{', re.MULTILINE)
RE_MAGIC_NUMBER = re.compile(r'(?<!\w)(?<!\.)(\d+\.\d{4,}|[1-9]\d{2,})(?!\.)(?!\w)') # More selective magic numbers
RE_IMPLEMENTATION_NAME = re.compile(r'\b\w+(?:_m|_km|_rad|_deg|_vec|_raw)\b') # Focus on variables
RE_RAW_DOUBLE = re.compile(r'\bdouble\b\s+(\w+)\s*[=;]')
RE_BOOLEAN_PARAMS = re.compile(r'\bbool\s+\w+')
RE_EMPTY_IF = re.compile(r'if\s*\(!?\w+\)\s*\{\s*\}')

def analyze_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    violations = []
    lines = content.splitlines()

    # Rule 1 & Rule 4 & Rule 8: Parser for functions
    # Note: Highly simplified C++ parsing
    stack = 0
    in_func = False
    func_start = 0
    func_name = ""
    
    # Pre-parse comments for Rule 3
    comments = []
    for i, line in enumerate(lines):
        match = re.search(r'//\s*(.*)', line)
        if match:
            comments.append((i, match.group(1).strip().lower()))

    # Find functions
    for i, line in enumerate(lines):
        # Very crude function detector
        if not in_func:
            match = re.search(r'(\w+)\s*\(([^)]*)\)\s*(?:const\s*)?\{', line)
            if match and not line.strip().startswith('//'):
                in_func = True
                func_start = i
                func_name = match.group(1)
                params = match.group(2)
                
                # Rule 8: Boolean Flags
                bool_count = len(RE_BOOLEAN_PARAMS.findall(params))
                if bool_count > 1:
                    violations.append((i + 1, "R8", f"Function '{func_name}' has {bool_count} boolean flags."))
                
                stack = line.count('{') - line.count('}')
                if stack == 0: # One-liner
                    in_func = False
        else:
            stack += line.count('{')
            stack -= line.count('}')
            if stack == 0:
                in_func = False
                length = i - func_start
                # Rule 1: Length
                if length > 25: # Allowing some headroom over 7-9 lines for C++
                    violations.append((func_start + 1, "R1", f"Function '{func_name}' is {length} lines long."))

        # Rule 5: Magic Numbers (Scan outside comments/strings)
        clean_line = re.sub(r'//.*', '', line)
        clean_line = re.sub(r'".*?"', '', clean_line)
        for m in RE_MAGIC_NUMBER.finditer(clean_line):
            violations.append((i + 1, "R5", f"Magic number '{m.group()}' detected in line."))

        # Rule 2: Implementation Names in variables (exclude method calls like .to_deg())
        for m in RE_IMPLEMENTATION_NAME.finditer(clean_line):
            name = m.group()
            # Check if it looks like a variable assignment or declaration, not a method call
            if not re.search(r'\.' + name + r'\(', clean_line):
                violations.append((i + 1, "R2", f"Potential implementation name '{name}' used for variable."))

        # Rule 6: Raw Double usage (exclude common indexes like i, j, k or results)
        for m in RE_RAW_DOUBLE.finditer(clean_line):
            var_name = m.group(1)
            if var_name not in ['i', 'j', 'k', 'dt', 'it', 'result', 'val']:
                violations.append((i + 1, "R6", f"Variable '{var_name}' uses raw double instead of physical type."))

        # Rule 7: Silent Fail
        if RE_EMPTY_IF.search(line):
            violations.append((i + 1, "R7", "Silent failure: empty IF block detected."))

    # Rule 3: Comments repeating code
    for i, comment_text in comments:
        if i + 1 < len(lines):
            next_line = lines[i+1].lower()
            # If comment words appear almost all in next line
            words = [w for w in comment_text.split() if len(w) > 3]
            if words and all(w in next_line for w in words):
                violations.append((i + 1, "R3", f"Comment repeats code: '{comment_text}'"))

    return violations

def main():
    if len(sys.argv) < 2:
        print("Usage: python code_fits_head_checker.py <dir_or_file>")
        return

    target = sys.argv[1]
    all_violations = {}

    if os.path.isfile(target):
        files = [target]
    else:
        files = []
        for root, _, fs in os.walk(target):
            for f in fs:
                if f.endswith(('.cpp', '.hpp', '.h')):
                    files.append(os.path.join(root, f))

    for f in files:
        results = analyze_file(f)
        if results:
            all_violations[f] = results

    # Print Summary
    print("\n" + "="*80)
    print(" CODE FITS IN YOUR HEAD - VIOLATION REPORT")
    print("="*80)
    
    total = 0
    for f, v_list in all_violations.items():
        rel_path = os.path.relpath(f)
        print(f"\nFILE: {rel_path}")
        for line, rule, msg in v_list:
            print(f"  [{rule}] Line {line:4}: {msg}")
            total += 1
            
    print("\n" + "="*80)
    print(f" TOTAL VIOLATIONS: {total}")
    print("="*80)

if __name__ == "__main__":
    main()
