package BackUp;

import java.util.ArrayList;
import java.util.List;

/**
 * cartesian
 *
 * @author ruiqi.zss@alibaba-inc.com
 * @date 2019/3/28
 */
public class CartesianSet<T> {

    private final T[][] source;

    private final long count;

    public CartesianSet(T[][] source) {
        this.source = source;
        long total = 1;
        for (T[] array : source) {
            total *= array.length;
        }
        this.count = total;
    }

    public long getCount() {
        return count;
    }

    public List<T> get(int index) {
        if (index < 0 || count <= index) { return null; }
        List<T> result = new ArrayList<T>(this.source.length);

        long weight = 1;
        for (T[] row : this.source) {
            long times = index / weight;
            int column = (int)(times % row.length);
            result.add(row[column]);
            weight *= row.length;
        }
        return result;
    }
}