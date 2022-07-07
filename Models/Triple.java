package Models;

public class Triple<U, V, Z>
{
    public final U first;       // the first field of a pair
    public final V second;      // the second field of a pair
    public final Z third;

    // Constructs a new pair with specified values
    public Triple(U first, V second, Z third)
    {
        this.first = first;
        this.second = second;
        this.third = third;
    }

    @Override
    // Checks specified object is "equal to" the current object or not
    public boolean equals(Object o)
    {
        if (this == o) {
            return true;
        }

        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        Triple<?, ?, ?> t = (Triple<?, ?, ?>) o;

        // call `equals()` method of the underlying objects
        if (!first.equals(t.first)) {
            return false;
        }
        if (!second.equals(t.second)) {
            return false;
        }
        return third.equals(t.third);
    }

    @Override
    // Computes hash code for an object to support hash tables
    public int hashCode()
    {
        // use hash codes of the underlying objects
        return 31 * first.hashCode() + second.hashCode() + third.hashCode();
    }

    @Override
    public String toString() {
        return "(" + first + ", " + second + ", " + third + ")";
    }

    // Factory method for creating a typed Pair immutable instance
    public static <U, V, Z> Triple <U, V, Z> of(U a, V b, Z c)
    {
        // calls private constructor
        return new Triple<>(a, b, c);
    }
}