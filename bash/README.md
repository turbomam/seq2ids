# `pgsql2sqlite`

## Script Usage

To run the script, use the following command:

```
poetry run sh bash/pgsql2sqlite.sh [username] [port_number] [table1,table2,table3,etc.] [dir_name]
```

## Arguments Description

`[username]`: This is the first argument, which as the name implies is the username to your SSH tunnel to LBL merlot.

`[port_number]`: This is the port number of the SSH tunnel on your local machine on which you are receving packets from LBL merlot.

`[table1,table2,table3,etc.]`: This is a list of comma separted table names that the script accepts. You will find output tables with corresponding names in the _dir_name_ folder.

`[dir_name]`: The directory where you want the script to dump SQLite exports. This is an optional argument. If you do not pass a value for it, it will dump all exports to a folder called `sqlite_dump`.
