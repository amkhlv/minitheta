Description
===========

This is a tool for manipulations with polinomials of fermionic variables.
It can add, multiply and differentiate polynomials, and solve linear equations. 

It was written to solve some linear algebra problems which arize in the compuation of
the [cohomology of the sum of two pure spinors](http://arxiv.org/abs/1301.3353)

The ScalaDoc is available [on my website](http://andreimikhailov.com/scaladoc/index.html)

I have not been working on it since; it runs on __Scala 2.9__

Running
=======

Creating jar file
-----------------

Go to that directory where you have cloned the minitheta.

    mkdir lib
    scalac -d lib `find -type f`
    cd lib
    jar cf minitheta.jar `find -name "*.class"`

Starting the REPL
-----------------

    JAVA_OPTS="-Xmx2g -Xss40m" scala -cp minitheta.jar

Now you are inside the Scala REPL. Which computation do you want to run? For example, you can say:

    com.andreimikhailov.brst.TLMTLT.main(Array())

Memory issues
-------------

Notice that we used the Java options `-Xmx2g` (this is heap size) and `-Xss40m` (this is the stack size).
But not every computer has 2Gb of memory, and on the other hand even 2Gb will not be enough for more complicated calculations.
So, this parameter has to be adjusted. 
