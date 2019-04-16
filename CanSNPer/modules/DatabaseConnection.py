import sys
import os
import sqlite3

class ConnectionError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class DatabaseConnection(object):
	"""docstring for DatabaseConnection"""
	def __init__(self, database, verbose=False):
		super().__init__()
		self.verbose = verbose
		self.database = database
		if not os.path.exists(self.database):
			if self.verbose:
				print("python modules/database/CreateDatabase.py {database}".format(database=self.database))
			os.system("python modules/database/CreateDatabase.py {database}".format(database=self.database))
		self.conn = self.connect(self.database)
		self.cursor = self.create_cursor(self.conn)

	def connect(self,database):
		'''Create database connection'''
		try:
			conn = sqlite3.connect(database)
			if self.verbose:
				print("{database} opened successfully.".format(database=database))

			return conn
		except Exception as e:
			sys.stderr.write(str(e))
		raise ConnectionError("Count not connect to the database {database} see above message for details!".format(database=database))
		return None

	def disconnect(self):
		self.conn.close()

	def create_cursor(self,conn):
		'''Create a db cursor'''
		try:
			cursor = conn.cursor()
			if self.verbose:
				print("cursor created.")
			return cursor
		except Exception as e:
			sys.stderr.write(str(e))
		return None

	def commit(self):
		self.conn.commit()

	def query(self,query,values = False, cursor=False):
		'''    The query function controls  the error handling of sqlite3 execute'''
		res = False
		if self.verbose: print(query)
		if not cursor:
			cursor = self.cursor
		try:
			if values:
				return cursor.execute(query,values)
			else:
				return cursor.execute(query)
		except Exception as e:
			if "UNIQUE constraint failed" not in str(e):
				## UNIQUE constraint is an accepted error as it keeps multiple edges from being added
				print("Error in DatabaseConnection query")
				if self.verbose: print(query)
				sys.stderr.write(str(e))
			return(e)

	def insert(self,data,table):
		'''Insert function
				data is a dictionary with keys matching
				table columns with respective value to be inserted
		'''
		INSERT_QUERY = '''
			INSERT INTO {table}({columns})
			VALUES ({values})
		'''

		columns = data.keys()
		values = tuple([data[col] for col in columns])
		insertStr = INSERT_QUERY.format(
				columns=",".join(columns),
				values=','.join(["?" for x in values]),
				table=table
		)
		return self.query(insertStr,insert_val=values)

class DatabaseFunctions(object):
	"""docstring for DatabaseFunctions."""
	def __init__(self, database, verbose=False):
		super().__init__()
		self.verbose = verbose
		if self.verbose: print("Load DatabaseFunctions")
		### store DatabaseConnection object reference
		self.database = database

	'''Set functions'''
	def set_database(self, database):
		'''Change the database object default in class'''
		self.database = database

	'''Get functions of class'''
	def get_taxid_base(self):
		'''Fetch the next incremental node from the current database'''
		QUERY = "SELECT MAX(id) AS max_id FROM nodes"
		return self.database.query(QUERY).fetchone()[0]+1

	def get_genomes(self, database=False,limit=0):
		'''Get the list of genomes in the database'''
		## This is a many to many relation, so all genomes has to be put in a set for each taxonomy id
		genomeDict = {}
		QUERY = '''SELECT id,genome FROM {table}'''.format(table="genomes")
		if limit > 0:
			QUERY += " LIMIT {limit}".format(limit=limit)
		print(QUERY)
		for id,genome in database.query(QUERY).fetchall():
			genomeDict[genome] = id
		return genomeDict

	def get_all(self, database=False, table=False):
		'''Get full table from table'''
		QUERY = '''SELECT * FROM {table}'''.format(table)
		return database.query(QUERY).fetchall()

	def get_nodes(self, database=False,col=False):
		'''Retrieve the whole node info table of the database to decrease the number of database calls!'''
		nodeDict = {}
		QUERY = '''SELECT id,name FROM nodes'''
		if not database:
			database = self.database
		for node in database.query(QUERY).fetchall():
			nodeDict[node[0]] = node[1]
			if col == 1:
				continue
			nodeDict[node[1]] = node[0]
		return nodeDict

	def get_links(self, nodes,database=False):
		'''This function returns all links in the given database'''
		QUERY = '''SELECT parent,child FROM tree WHERE parent in ({nodes}) OR child in ({nodes})'''.format(nodes=",".join(map(str,nodes)))
		if not database:
			database = self.database
		links = database.query(QUERY).fetchall()
		return links

	'''Add functions of class'''
	def add_node(self, description, id=False):
		'''Add node to tree'''
		info = { "name": description }
		### If ID is supplied skip autoincrement and add specific ID
		if id:
			info["id"] = id
		taxid_base = self.database.insert(info, table="nodes")
		return taxid_base

	def add_link(self, child, parent,table="tree"):
		'''Add relationship in tree'''
		info = {
			"child": child,
			"parent": parent
		}
		return self.database.insert(info, table="tree")

	def add_genome(self, id, genome):
		'''Add genome annotation to nodes'''
		info = {
			"id": id,
			"genome": genome
		}
		return self.database.insert(info, table="genomes")

	def add_links(self,links, table="tree",hold=False):
		'''Add links from a list to tree'''
		added_links = 0
		for parent,child in links:
			res = self.add_link(parent,child,table=table)
			### Check if the link already exist in the database, this overlap may occur when a large new branch is added
			if "UNIQUE constraint failed" not in str(res):
				added_links +=1
		## Commit changes
		if not hold:
			self.database.commit()
		return added_links

	def add_nodes(self,nodes, table="tree",hold=False):
		'''Add nodes from a list of nodes'''
		added_nodes = 0
		for node in nodes:
			res = self.add_node(node)
			### Check if the link already exist in the database, this overlap may occur when a large new branch is added
			if "UNIQUE constraint failed" not in str(res):
				added_nodes +=1
		## Commit changes
		if not hold:
			self.database.commit()
		return added_nodes

	'''Delete functions of class'''
	def delete_links(self,links, table="tree",hold=False):
		'''This function deletes all links given in links'''
		QUERY = "DELETE FROM {table} WHERE parent = {parent} AND child = {child}"
		for parent,child in links:
			res = self.database.query(QUERY.format(table=table, parent=parent, child=child))
		## Commit changes
		if not hold:
			self.database.commit()

	def delete_nodes(self, nodes, table="nodes",hold=False):
		'''This function deletes all nodes given in nodes'''
		QUERY = "DELETE FROM {table} WHERE id = {node}"
		for node in nodes:
			res = self.database.query(QUERY.format(table=table, node=node))
		## Commit changes
		if not hold:
			self.database.commit()

	def num_rows(self,table):
		'''Return the number of rows in a table'''
		QUERY = '''SELECT Count(*) FROM {table}'''
		return self.database.query(QUERY.format(table=table)).fetchall()[0][0]
