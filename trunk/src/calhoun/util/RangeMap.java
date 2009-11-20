/*
 * The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This software and its documentation are copyright 2003 by
 * the Broad Institute/Massachusetts Institute of Technology. All rights are reserved. This software is supplied
 * without any warranty or guaranteed support whatsoever. Neither the Broad Institute nor MIT can be responsible for
 * its use, misuse, or functionality.
 */
package calhoun.util;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * A rangeMap contains a set of intervals and maps each interval to a set of objects that exist in that interval. It
 * allows very fast lookups of range queries.
 *
 * An interval is closed, that is, the interval from 20-30 includes both 20 and 30.
 */
public class RangeMap implements Serializable {

	public static final long serialVersionUID = 413339879647819935L;
	
	/**
	 * Imagine a number line. Take the set of start and stop coordinates for each object in the RangeMap and place them
	 * on the number line. This divides the number line into a series of intervals. Each interval can be associated
	 * with a distinct list of objects that contain all points in that interval. The RangeMap is a sorted map that keys
	 * the start of the interval to this distinct list of objects.
	 * <p>
	 * Objects that fully contain the points in more than one interval will be contained in multiple lists. One list
	 * for each interval they contain.
	 */
	SortedMap map = new TreeMap();

	/**
	 * Keeps track of the start and end for each object in the map so we can
	 * remove them easily if necessary.
	 *
	 * Maps Object -> int[2]. the first element in the array is the start of
	 * the key and the second element is the end of the key.
	 */
	Map objectList = new HashMap();

	public RangeMap() {
		super();
	}

	public int size() {
		return objectList.size();
	}

	public Set values() {
		return objectList.keySet();
	}

	/** Returns the start of the first non-empty interval.  If the RangeMap has no entries, returns 0. */
	public int getStart() {
		return size() == 0 ? 0 : ((Integer) map.firstKey()).intValue();
	}
	
	/** Returns the end of the last interval */
	public int getStop() {
		return size() == 0 ? 0 : ((Integer) map.lastKey()).intValue() - 1;
	}

	/**
	 * An interval is a convenient data structure used to return a section of a RangeMap. It contains a start, a stop,
	 * and a set of elements that occur in that range.
	 */
	public static class Interval {
		public int start;
		public int stop;
		public Set elements;
		
		public Interval() {}
		
		public Interval(int start, int stop, Set elements)
		{
			this.start = start;
			this.stop = stop;
			this.elements = new HashSet(elements);
		}
		
		@Override
		public boolean equals(Object obj)
		{
			if (!(obj instanceof Interval)) {
				return false;
			}

			Interval interval = (Interval) obj;

			if (this.start != interval.start) {
				return false;
			} else if (this.stop != interval.stop) {
				return false;
			} else if (!this.elements.containsAll(interval.elements)) {
				return false;
			} else if (!interval.elements.containsAll(this.elements)) {
				return false;
			}
		
			return true;
		}

		@Override
		public String toString() {
			String result = "";
			for (Iterator iter = elements.iterator(); iter.hasNext(); ){
				result += iter.next().toString();
			}
			return "interval " + start +" -> " +  stop + " : "+ result + "\n";
		}
		
		public int getLength()
		{
			return 1 + stop - start;
		}
	}

	public SortedMap getMap()
	{
		return Collections.unmodifiableSortedMap(map);
	}
	
	public Map getObjectList()
	{
		return Collections.unmodifiableMap(objectList);
	}
	
	public List getDisjointRegions()
	{
		return getRegions(true);
	}
	
	public List getRegions()
	{
		return getRegions(false);
	}
	
	/**
	 * Returns a list of Interval objects that contain overlapping elements.
	 *
	 * @param splitDisjointRegions if true, return disjoint regions,
	 * if false, return "old-style" regions. Consider the following:
	 *
	 *    <---------Feature A------->
	 *                               <----------Feature B---------->
	 *    <---------Feature C------->
	 *
	 * - If we group these into a single region, then any two groups
	 *   returned by this method will be separated by at least one gap
	 *   (ie an area that contains no elements at all). This is more
	 *   convenient.
	 * - If we group them into two regions, then every element in
	 *   a region will overlap at least one other element in that
	 *   region (if the region contains more than one element).
	 *   This is more mathematically precise.
	 *
	 * There are good reasons to want either behavior so we support both.
	 * Disjoint regions are more mathematically consistent but may lie
	 * immediately adjacent to each other. Nondisjoint regions are better
	 * separated, but "old-style" regions may contain an element that
	 * does not overlap anything else within it.
	 *
	 * No element in a region will ever overlap an element outside the region.
	 *
	 * This method returns a List of Interval objects. The Interval objects
	 * returned will be in increasing order of their placement on the RangeMap.
	 *
	 * @see #getIntervals()
	 */
	public List getRegions(boolean splitDisjointRegions) {
		List list = new ArrayList();
		Interval currentInterval = null;
		// Iterate through the map (in order). As long as there are intervals with at least one entry, add that entry
		// to the current cluster.
		// When you hit an interval with no values, finish the interval to return.
		Iterator it = map.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry entry = (Map.Entry) it.next();
			List currentValue = (List) entry.getValue();
			int currentKey = ((Integer) entry.getKey()).intValue();
			if (currentValue.size() == 0) {
				Assert.a(currentInterval != null,
				  "Multiple 0 intervals in a row were detected.");
				currentInterval.stop = currentKey - 1;
				list.add(currentInterval);
				currentInterval = null;
			} else {
				if (currentInterval != null) {
					/*
					 * If two intervals are adjacent but disjoint (no element
					 * spans them) then return them as two groups, not one.
					 */
					boolean disjoint = true;
					for (Iterator i = currentValue.iterator(); i.hasNext(); ) {
						if (currentInterval.elements.contains(i.next())) {
							disjoint = false;
						}
					}
					if (disjoint && splitDisjointRegions) {
						currentInterval.stop = currentKey - 1;
						list.add(currentInterval);
						currentInterval = null;
					}
				}
			
				if (currentInterval == null) {
					currentInterval = new Interval();
					currentInterval.start = currentKey;
					currentInterval.elements = new HashSet();
				}
				currentInterval.elements.addAll(currentValue);
			}
		}
		Assert.a(currentInterval == null, "RangeMap did not end with a 0 entry.");
		return list;
	}

	/**
	 * Returns the individual intervals in the rangeMap.
	 *
	 * An interval is an area of a range map where every point in that
	 * area is overlapped by the same set of elements. Consider this:
	 *
	 *    <-----Feature A----->           <---Feature C--->
	 *            <---------Feature B--------->
	 *
	 *    |   1   |     2     |     3     | 4 |     5     |
	 *
	 * This would be represented by one region but counts as five
	 * intervals: interval 1 contains A only, interval 2 contains A & B,
	 * interval 3 holds B only, and so on.
	 *
	 * Callers should be prepared for Intervals that contain no elements.
	 *
	 * This method returns a List of Interval objects. The Interval objects
	 * returned will be in increasing order of their placement on the RangeMap.
	 *
	 * @see #getRegions()
	 */
	public List getIntervals()
	{
		List list = new ArrayList();
		Interval currentInterval = null;
		Iterator it = map.entrySet().iterator();

		while (it.hasNext()) {
			Map.Entry entry = (Map.Entry) it.next();
			List currentValue = (List) entry.getValue();
			int currentKey = ((Integer) entry.getKey()).intValue();
			
			if (currentInterval != null) {
				currentInterval.stop = currentKey - 1;
				list.add(currentInterval);
			}
			
			currentInterval = new Interval();
			currentInterval.start = currentKey;
			currentInterval.elements = new HashSet();
			if (currentValue != null) {
				currentInterval.elements.addAll(currentValue);
			}
		}
		
		Assert.a(currentInterval.elements.size() == 0,
		  "RangeMap did not end with an empty entry");
		return list;
	}

	/** Like getIntervals, but allows you to define a start and stop.  First and last interval are guaranteed to start and stop at the endpoints. */
	public List getIntervals(int start, int stop)
	{
		List list = new ArrayList();
		Interval currentInterval = null;
		Iterator it = getSubMap(start, stop+1).entrySet().iterator();

		while (it.hasNext()) {
			Map.Entry entry = (Map.Entry) it.next();
			List currentValue = (List) entry.getValue();
			int currentKey = ((Integer) entry.getKey()).intValue();
			
			if(currentInterval == null && start < currentKey) {
				// We started before the beginning of the range map, add in a dummy empty interval
				currentInterval = new Interval();
				currentInterval.start = start;
				currentInterval.elements = Collections.EMPTY_SET;
			}

			if (currentInterval != null) {
				currentInterval.stop = currentKey - 1;
				list.add(currentInterval);
			}

			currentInterval = new Interval();
			currentInterval.start = currentKey > start ? currentKey : start;
			if (currentValue == null) {
				currentInterval.elements = Collections.EMPTY_SET;
			}
			else {
				currentInterval.elements = new HashSet(currentValue);
			}
		}

		if(currentInterval != null) {
			currentInterval.stop = stop;
			list.add(currentInterval);
		}
		return list;
	}	
			
	/** Returns the empty intervals.  A start and stop are required to bound the edges of the first and last region.  These bounds can be within the existin region.
	 * The returned intervals are closed.
	*/
	public List getEmptyIntervals(int start, int totalSize) {
		Assert.a(start <= totalSize, "Start (",new Integer(start),") must be less than stop (",new Integer(totalSize),")");
		List list = new ArrayList();
		Interval currentInterval = null;
		
		boolean beforeStart = true;
		Iterator it = map.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry entry = (Map.Entry) it.next();
			List currentValue = (List) entry.getValue();
			int currentKey = ((Integer) entry.getKey()).intValue();
			if(currentKey > totalSize)
				break;
			if (currentValue.size() == 0) {
				Assert.a(beforeStart == false, "Map started with a 0 entry.");
				Assert.a(currentInterval == null, "Multiple 0 intervals in a row were detected.");
				if(currentKey < start) {
					beforeStart = true;
				}
				else {
					currentInterval = new Interval();
					currentInterval.start = currentKey;
				}
			} else {
				if(beforeStart == true) {
					beforeStart = false;
					if(start < currentKey) {
						currentInterval = new Interval();
						currentInterval.start = start;
					}
				}
				if(currentInterval != null) {
					currentInterval.stop = currentKey - 1;
					list.add(currentInterval);
					currentInterval = null;
				}
			}
		}
		if(beforeStart == true) {
			Assert.a(currentInterval == null);
			currentInterval = new Interval();
			currentInterval.start = start;
		}
		if(currentInterval != null) {
			currentInterval.stop = totalSize;
			list.add(currentInterval);
		}
		return list;
	}

	/**
	 * Returns objects in the range map that overlap this range.
	 *
	 * query:      *-----------------*      returned
	 *          *------------*              returned
	 *                *-------------*       returned
	 *           *---------------------*    returned
	 *        *----*                        returned
	 *                               *-*    returned
	 */
	public Set find(int argLow, int argHigh) {
		Iterator it = getSubMap(argLow, argHigh+1).values().iterator();
		Set ret = new HashSet();
		while (it.hasNext()) {
			ArrayList l = (ArrayList) it.next();
			ret.addAll(l);
		}
		return ret;
	}

	/**
	 * Returns objects in the range map that are fully contained by the range.
	 *
	 * query:      *-----------------*      returned
	 *          *------------*              not returned
	 *                *-------------*       returned
	 *           *---------------------*    not returned
	 *        *----*                        not returned
	 *                               *-*    not returned
	 */
	public Set findContained(int argLow, int argHigh) {
		Iterator it = find(argLow, argHigh).iterator();

		Set ret = new HashSet();
		while (it.hasNext()) {
			Object o = it.next();
			Integer [] bounds = (Integer[]) objectList.get(o);
			if(bounds[0].intValue() >= argLow && bounds[1].intValue() <= argHigh+1) {
				ret.add(o);
			}
		}
		return ret;
	}
	
	/**
	 * Returns objects in the range map that fully contain the range.
	 *
	 * query:      *-----------------*      returned
	 *          *------------*              not returned
	 *                *-------------*       not returned
	 *           *---------------------*    returned
	 *        *----*                        not returned
	 *                               *-*    not returned
	 */
	public Set findContaining(int argLow, int argHigh) {
		Iterator it = find(argLow, argHigh).iterator();

		Set ret = new HashSet();
		while (it.hasNext()) {
			Object o = it.next();
			Integer [] bounds = (Integer[]) objectList.get(o);
			if(bounds[0].intValue() <= argLow && bounds[1].intValue() >= argHigh+1) {
				ret.add(o);
			}
		}
		return ret;
	}
	
	/** Returns true if the map has an entry in the specfied range. */
	public boolean hasEntry(int argLow, int argHigh) {
		Iterator it = getSubMap(argLow, argHigh+1).values().iterator();

		// If there is no entry or a single empty entry return false.  Otherwise, return true
		if(!it.hasNext() || ((ArrayList) it.next()).size() == 0 && !it.hasNext()) {
			return false;
		}
		else
			return true;
	}
	
	/** Returns the submap from argLow+1 included to argHigh excluded. */
	protected SortedMap getSubMap(int argLow, int argHigh) {
		if (argLow > argHigh)
			throw new IllegalArgumentException(
				"Low end of range (" + argLow + ") is greater than high end (" + argHigh + ")");
		Integer lowBound = new Integer(argLow + 1);
		Integer high = new Integer(argHigh);

		SortedMap lowMap = map.headMap(lowBound);
		if (lowMap.size() != 0)
			lowBound = (Integer) lowMap.lastKey();
		return map.subMap(lowBound, high);
	}

	/**
	 * Adds the object o into the map with interval lo - high. The interval is inclusive, with both the high and low being part of the interval.
	 */
	public void add(int argLow, int argHigh, Object o) {
		if (argLow > argHigh)
			throw new IllegalArgumentException(
				"Low end of range (" + argLow + ") is greater than high end (" + argHigh + ")");
		if (objectList.containsKey(o)) {
			throw new IllegalArgumentException("RangeMap already contains " + o);
		}
		Integer low = new Integer(argLow);
		Integer high = new Integer(argHigh+1);
		objectList.put(o, new Integer[] { low, high });

		makeSplitAt(low);
		makeSplitAt(high);

		// Take the new submap and add o to every element
		Iterator it = map.subMap(low, high).values().iterator();
		while (it.hasNext()) {
			ArrayList l = (ArrayList) it.next();
			l.add(o);
		}
	}
	

	/** Returns true if the object is in the RangeMap */
	public boolean contains(Object o) {
		return objectList.get(o) != null;
	}

	/**
	 * Removes the object o from the map.
	 */
	public void remove(Object o) {
		Integer[] bounds = (Integer[]) objectList.get(o);
		if (bounds == null) {
			throw new IllegalArgumentException("RangeMap does not contain: " + o);
		}
		objectList.remove(o);
		// Remove the object from each interval it contained.
		Iterator it = map.subMap(bounds[0], bounds[1]).values().iterator();
		while (it.hasNext()) {
			ArrayList l = (ArrayList) it.next();
			l.remove(o);
		}
		/*
		 * Check to see if we can remove an interval. This occurs if there is no change from the
		 */
		mergeAt(bounds[0]);
		mergeAt(bounds[1]);
	}

	private void mergeAt(Integer bound) {
		List startList = (List) map.get(bound);
		SortedMap lowMap = map.headMap(bound);
		if (lowMap.size() == 0) {
			if (startList.size() == 0)
				map.remove(bound);
		} else {
			List previousList = (List) lowMap.get(lowMap.lastKey());
			if (previousList.equals(startList)) {
				map.remove(bound);
			}
		}
	}

	private void makeSplitAt(Integer bound) {
		// If there is already a split, we are done
		Object entry = map.get(bound);
		if (entry == null) {
			// Find the previous interval and duplicate it
			SortedMap headMap = map.headMap(bound);
			if (headMap.size() == 0) {
				map.put(bound, new ArrayList());
			} else {
				map.put(bound, new ArrayList((List) map.get(headMap.lastKey())));
			}
		}
	}

	@Override
	public String toString() {
		StringBuffer ret = new StringBuffer("Range Map (x-y:count) ");
		Integer lastKey = null;
		int lastCount = 0;
		Iterator it = map.keySet().iterator();
		while (it.hasNext()) {
			Integer key = (Integer) it.next();
			if (lastKey != null) {
				ret.append(lastKey + "-" + (key.intValue() - 1) + ":" + lastCount + " ");
			}
			lastKey = key;
			lastCount = ((List) map.get(key)).size();
		}
		ret.append(lastKey + "-end:" + lastCount);
		return ret.toString();
	}
}
