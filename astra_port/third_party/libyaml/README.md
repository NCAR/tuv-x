# Bundled LibYAML

Unmodified C sources and headers from [LibYAML 0.2.5](https://github.com/yaml/libyaml/tree/0.2.5),
commit `2c891fc7a770e8ba2fec34fc6b545c672beb37e6`, under the included MIT license.
`provenance.json` records SHA-256 hashes of every copied upstream file.

The TurboTUV build includes `namespace.h` before these sources to prefix all
LibYAML identifiers with `tt_`. Its public types do not appear in the TurboTUV
API. This avoids collisions with host applications using another LibYAML copy.
The bundled parser requires no runtime library, package manager or download.
