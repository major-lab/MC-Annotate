#include "StemConnection.h"

namespace annotate 
{
	StemConnection::StemConnection() 
	{
		mpStem = NULL; 
		meConnection = Stem::eUNDEFINED_CONNECTION;
	}
		
	StemConnection::StemConnection(
		const Stem& aStem, 
		const Stem::enConnection& aeConnect)
	{
		mpStem = &aStem;
		meConnection = aeConnect;
	}
	
	int StemConnection::getDirection() const
	{
		int iDirection = 0;
		if(isValid())
		{
			switch(meConnection)
			{
			case Stem::eFIRST_STRAND_FRONT_PAIR:
				iDirection = -1;
				break;
			case Stem::eFIRST_STRAND_BACK_PAIR:
				iDirection = 1;
				break;
			case Stem::eSECOND_STRAND_FRONT_PAIR:
				if(Stem::eANTIPARALLEL == mpStem->getOrientation())
				{
					iDirection = 1;
				}
				else
				{
					iDirection = -1;
				}
				break;
			case Stem::eSECOND_STRAND_BACK_PAIR:
				if(Stem::eANTIPARALLEL == mpStem->getOrientation())
				{
					iDirection = -1;
				}
				else
				{
					iDirection = 1;
				}
				break;
			default:
				// TODO : Exception
				iDirection = 0;
			}
		}
		return iDirection;
	}
	
	mccore::ResId StemConnection::nextId() const
		throw(mccore::NoSuchElementException)
	{
		mccore::ResId id;
		if(NULL != mpStem && meConnection != Stem::eUNDEFINED_CONNECTION)
		{
			switch(meConnection)
			{
			case Stem::eFIRST_STRAND_FRONT_PAIR:
				id = mpStem->basePairs().front().rResId;
				break;
			case Stem::eSECOND_STRAND_FRONT_PAIR:
				id = mpStem->basePairs().front().fResId;
				break;
			case Stem::eFIRST_STRAND_BACK_PAIR:
				id = mpStem->basePairs().back().rResId;
				break;
			case Stem::eSECOND_STRAND_BACK_PAIR:
				id = mpStem->basePairs().back().fResId;
				break;
			default:
				std::string strMsg("Invalid connections have no nextId");
				throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
				break;
			}
		}
		else
		{
			std::string strMsg("Invalid connections have no nextId");
			throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
		}
		
		return id;
	}
	
	BasePair StemConnection::getPair() const
		throw(mccore::NoSuchElementException)
	{
		const BasePair* pBasePair = NULL;
		if(NULL != mpStem && meConnection != Stem::eUNDEFINED_CONNECTION)
		{
			switch(meConnection)
			{
			case Stem::eFIRST_STRAND_FRONT_PAIR:
			case Stem::eSECOND_STRAND_FRONT_PAIR:
				pBasePair = &mpStem->basePairs().front();
				break;
			case Stem::eFIRST_STRAND_BACK_PAIR:
			case Stem::eSECOND_STRAND_BACK_PAIR:
				pBasePair = &mpStem->basePairs().back();
				break;
			default:
				pBasePair = NULL;
			}
		}
		
		if(NULL == pBasePair)
		{
			std::string strMsg("Invalid connections have no connecting pairs");
			throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
		}
		return *pBasePair;		
	}
	
	bool StemConnection::connects(const StemConnection& aConnection) const
	{
		bool bConnects = false;
		if(isValid() && aConnection.isValid())
		{
			if(mpStem == aConnection.mpStem)
			{
				bConnects = (nextId() == aConnection.getResidue());
			}
		}
		return bConnects;
	}
	
	mccore::ResId StemConnection::getResidue() const 
	{ 
		return mpStem->getResidue(meConnection); 
	}
	
	bool StemConnection::isValid() const 
	{ 
		return (NULL != mpStem && meConnection != Stem::eUNDEFINED_CONNECTION);
	}
	
	bool StemConnection::operator== (const StemConnection &other) const 
	{
		bool bEqual = (mpStem == other.mpStem) 
			&& (meConnection == other.meConnection);
		return bEqual;
	}
};

