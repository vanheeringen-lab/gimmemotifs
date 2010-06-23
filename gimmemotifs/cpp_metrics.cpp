//  Copyright Joel de Guzman 2002-2004. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//  Hello World Example from the tutorial
//  [Joel de Guzman 10/9/2002]

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>
using std::cout;

char const* greet()
{
   return "hello, world";
}

char const* ten()
{
	boost::math::chi_squared mydist(3);
	double p = boost::math::cdf(mydist,.660);
	cout << "p-value " << p;
	return "10";
}

BOOST_PYTHON_MODULE(hello_ext)
{
    using namespace boost::python;
    def("greet", greet);
    def("ten", ten);
}

