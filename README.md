First the groovy migration script `migrate_annotations1to2.groovy`

There was a mix of transcript and mRNA after migration, so I changed every transcript to mRNA in the db:

```
update feature set class='org.bbop.apollo.MRNA' where class = 'org.bbop.apollo.Transcript';
```

Then the script to transfer attributes:

```
python migrate.py
```

Manually created the users on the new Apollo2 server with arrow (=python-apollo; we use ldap auth on our nginx proxy, so creating REMOTE_USER):

```
randomPass=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1)
DEFAULT_GROUP='annotators'

arrow users create_user --role user --metadata '{"INTERNAL_PASSWORD":"'$randomPass'"}' someone@xxxxxxx REMOTE_USER someone@xxxxxxx $randomPass
arrow users add_to_group $DEFAULT_GROUP someone@xxxxxxx
```

And finally the script to fix owners:

```
python owners.py | grep update > owners.sql && psql -h db -Upostgres < owners.sql
```

All these scripts use bcbio-gff and python-apollo modules: `pip install bcbio-gff git+https://github.com/galaxy-genome-annotation/python-apollo.git`
