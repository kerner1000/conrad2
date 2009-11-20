package calhoun.analysis.crf.executables;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import calhoun.util.FileUtil;

public class WeightCompare {
	public static void main(String[] args) throws Exception {
		String s = FileUtil.readFile("T:/crf/cnDT/data/train_50_5/name.dat");
		Set a = new HashSet();
		a.addAll(Arrays.asList(s.split("\n")));
		overlap(a, "T:/crf/cnDT/data/train_100_7/name.dat");
		overlap(a, "T:/crf/cnDT/data/train_1000_9/name.dat");
	}

	static void overlap(Set a, String name) throws Exception {
		String s = FileUtil.readFile(name);
		List<String> l = Arrays.asList(s.split("\n"));
		System.out.println(a.contains(a.iterator().next()));
		
		int i=0;
		for(String q : l) {
			if(a.contains(q)) {
				++i;
				System.out.println(q);
			}
		}
		System.out.println(i);
	}
}
