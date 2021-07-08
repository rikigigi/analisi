#ifndef MACROS_H
#define MACROS_H

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)

#endif // MACROS_H
