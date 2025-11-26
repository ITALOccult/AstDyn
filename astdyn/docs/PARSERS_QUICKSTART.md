# Parser Configurabili - Quick Start

## Uso Base

```cpp
#include "astdyn/io/IParser.hpp"

using namespace astdyn::io;

// Auto-detection del formato
auto parser = ParserFactory::create_orbit_parser("203.eq1");
auto elements = parser->parse("203.eq1");

// Accedi agli elementi
std::cout << "a = " << elements.semi_major_axis << " AU\n";
std::cout << "e = " << elements.eccentricity << "\n";
```

## Formati Supportati

| Tipo | Estensione | Auto-detect | Manuale |
|------|-----------|-------------|---------|
| OrbFit EQ1 | `.eq1` | ✅ | `"eq1"` |
| OrbFit RWO | `.rwo` | ✅ | `"rwo"` |

## Esempi

### 1. Parse con auto-detection
```cpp
auto parser = ParserFactory::create_orbit_parser("data.eq1");
auto elem = parser->parse("data.eq1");
```

### 2. Parse con formato esplicito
```cpp
auto parser = ParserFactory::create_orbit_parser("data.txt", "eq1");
auto elem = parser->parse("data.txt");
```

### 3. Check capabilities
```cpp
std::cout << parser->name() << "\n";
std::cout << parser->can_handle("test.eq1") << "\n";
```

## Documentazione Completa

Vedi [CONFIGURABLE_PARSERS.md](CONFIGURABLE_PARSERS.md) per dettagli su:
- Architettura del sistema
- Come aggiungere nuovi formati
- Integrazione con AstDyn Engine
- Note tecniche sulla conversione equinoziale→kepleriana
