print("""#define DECL_ARG(a,b) a b
#define NAME_ARG(a,b) b""")


for i in range(1,10):
    print ("#define DECL_CALL_BASE_{}(ret, fname".format(i),end="")
    for j in range(1,i+1):
        print (", arg{}".format(j),end="")
    print (") ret fname ( DECL_ARG arg1",end="")
    for j in range(2,i+1):
        print (", DECL_ARG arg{}".format(j),end="")
    print (") { return static_cast<T*>(this)-> fname ( NAME_ARG arg1",end="")
    for j in range(2,i+1):
        print (", NAME_ARG arg{}".format(j),end="")
    print(" ); }")
