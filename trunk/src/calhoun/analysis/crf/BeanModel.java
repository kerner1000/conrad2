package calhoun.analysis.crf;

import java.io.Serializable;
import java.util.List;

import calhoun.util.DenseBooleanMatrix2D;

/** a useful base class for creating model beans.  Defines internal Node and Edge class that can be used to define a hidden state model for the CRF.
 * The actual state diagram could be implemented through a derived class or through setup of the model in the XML model file.  Because XML configuration of
 * a state model is verbose and because the state diagram is central to most models, we recommend using a derived class to set up the model, as we have
 * done with gene calling.
 */
public class BeanModel extends CompositeFeatureManager implements ModelManager {
	private static final long serialVersionUID = -6641426879906871691L;

	/** a hidden state.  Each state has a name and states are numbered beginning from 0.
	 * The state index is used when the engine predicts hidden sequences.
	 */
	public static class Node implements Serializable {
		private static final long serialVersionUID = 1550966522391566799L;

		int index;
		String name;

		public Node() { }
		public Node(int index, String name) { 
			this.index = index;
			this.name = name;
		}
		
		public int getIndex() { return index; }
		public void setIndex(int index) { this.index = index; }
		public String getName() { return name; }
		public void setName(String name) { this.name = name; }
	}
	
	/** an edge in the hidden state diagram.  Presence of an edge indicates that a transition between two states is valid. *
	 */
	public static class Edge implements Serializable {
		private static final long serialVersionUID = 5429589777288279645L;

		Node from;
		Node to;

		public Edge() { }
		public Edge(Node from, Node to) { 
			this.from = from;
			this.to = to;
		}

		public Node getFrom() { return from; }
		public void setFrom(Node from) { this.from = from; }
		public Node getTo() { return to; }
		public void setTo(Node to) { this.to = to; }
	}
	
	protected List<Node> nodes;
	protected List<Edge> edges;

	/** sets the list of hidden states in this model.  This will usually be called by a derived class 
	 * or by the Spring configuration.
	 * @param nodes the list of hidden states in the model
	 */
	public void setNodes(List<Node> nodes) {
		this.nodes = nodes;
	}
	
	/** sets the list of legal state transitions in this model.  This will usually be called by a derived class 
	 * or by the Spring configuration.
	 * @param edges the list of valid state transitions in the model
	 */
	public void setEdges(List<Edge> edges) {
		this.edges= edges;
	}
	
	public int getNumStates() {
		return nodes.size();
	}
	
	public String getStateName(int state) {
		return nodes.get(state).name;
	}
	
	public int getStateIndex(String name) {
		for(Node n : nodes) {
			if(n.name.equals(name))
				return n.index;
		}
		throw new IllegalArgumentException();
	}

	public DenseBooleanMatrix2D getLegalTransitions() {
		int nStates = getNumStates();
		DenseBooleanMatrix2D ret = new DenseBooleanMatrix2D(nStates, nStates);
		for(Edge e : edges) {
			ret.setQuick(e.from.index, e.to.index, true);
		}
		return ret;
	}
}
