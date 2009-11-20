package calhoun.util;

/** Assertion class.  Like Java assert except it can't be shut off.
 * Assert.a(condition) or Assert.a(condition, msg) 
 * Some take multiple objects.  The advantage is that the string is only concatenated if the assert is fired.
*/
public final class Assert {
	private Assert() { }

	public static final void a(final boolean bool) {
		if(!bool)
			throw new CheckException();
	}

	public static final void a(final boolean bool, final String message) {
		if(!bool)
			throw new CheckException(message);
	}

	public static final void a(final boolean bool, final Object obj1, final Object obj2) {
		if(!bool)
			throw new CheckException("" + obj1 + obj2);
	}

	public static final void a(final boolean bool, final Object obj1, final Object obj2, final Object obj3) {
		if(!bool)
			throw new CheckException("" + obj1 + obj2 + obj3);
	}

	public static final void a(final boolean bool, final Object obj1, final Object obj2, final Object obj3, final Object obj4) {
		if(!bool)
			throw new CheckException("" + obj1 + obj2 + obj3 + obj4);
	}

	public static final void a(final boolean bool, final Object obj1, final Object obj2, final Object obj3, final Object obj4, final Object obj5) {
		if(!bool)
			throw new CheckException("" + obj1 + obj2 + obj3 + obj4 + obj5);
	}

	public static final void a(final boolean bool, final Object obj1, final Object obj2, final Object obj3, final Object obj4, final Object obj5, final Object obj6) {
		if(!bool)
			throw new CheckException("" + obj1 + obj2 + obj3 + obj4 + obj5 + obj6);
	}

	public static final void a(final boolean bool, final Object obj1, final Object obj2, final Object obj3, final Object obj4, final Object obj5, final Object obj6, final Object obj7) {
		if(!bool)
			throw new CheckException("" + obj1 + obj2 + obj3 + obj4 + obj5 + obj6 + obj7);
	}

	public static final void a(final boolean bool, final Object obj1, final Object obj2, final Object obj3, final Object obj4, final Object obj5, final Object obj6, final Object obj7, final Object obj8) {
		if(!bool)
			throw new CheckException("" + obj1 + obj2 + obj3 + obj4 + obj5 + obj6 + obj7 + obj8);
	}

	public static final void a(final boolean bool, final Object obj1, final Object obj2, final Object obj3, final Object obj4, final Object obj5, final Object obj6, final Object obj7, final Object obj8, final Object obj9) {
		if(!bool)
			throw new CheckException("" + obj1 + obj2 + obj3 + obj4 + obj5 + obj6 + obj7 + obj8 + obj9);
	}

	public static final void a(final boolean bool, final Object obj1, final Object obj2, final Object obj3, final Object obj4, final Object obj5, final Object obj6, final Object obj7, final Object obj8, final Object obj9, final Object obj10) {
		if(!bool)
			throw new CheckException("" + obj1 + obj2 + obj3 + obj4 + obj5 + obj6 + obj7 + obj8 + obj9 + obj10);
	}
}