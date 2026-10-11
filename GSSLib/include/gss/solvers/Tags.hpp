#pragma once

namespace gss {

struct PicardTag {};
struct NewtonTag {};
// EXTENSION POINT: add new Tag specialization here
// See docs/extending.rst for instructions

struct CholeskySparseTag {};
struct AMGTag {};
// EXTENSION POINT: add new Tag specialization here
// See docs/extending.rst for instructions

struct ZeroGuessTag {};
struct FadeevPerturbedTag {};
// EXTENSION POINT: add new Tag specialization here
// See docs/extending.rst for instructions

struct UniformVoronoiTag {};
// EXTENSION POINT: add new Tag specialization here
// See docs/extending.rst for instructions

} // namespace gss
