#pragma once

#define offset_of(st, m) ((size_t)((char *)&((st *)0)->m - (char *)0))
